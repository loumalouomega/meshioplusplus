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

// The frame-to-text and frame-to-file half of operations/render.hpp: the
// terminal cell encodings, the graphics protocols, PNG, asciicast and
// write_snapshot.

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/render.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"

// Project includes (private, not installed)
#include "../detail/png_write.hpp"
#include "../detail/raster.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kTxPrefix = "meshio++: render: ";

using TxRgb = std::array<int, 3>;

// ---------------------------------------------------------------------------
// Cell grids and glyphs
// ---------------------------------------------------------------------------

bool tx_is_cell(TextEncoding Encoding) {
    return Encoding == TextEncoding::HalfBlock || Encoding == TextEncoding::Quadrant ||
           Encoding == TextEncoding::Sextant || Encoding == TextEncoding::Braille ||
           Encoding == TextEncoding::Ascii;
}

// Subpixels per cell, columns by rows. Subpixel index is row-major, which is
// also the numbering of the sextant code points (1 top-left, 2 top-right, ...).
std::array<int, 2> tx_grid(TextEncoding Encoding) {
    switch (Encoding) {
        case TextEncoding::HalfBlock:
            return {1, 2};
        case TextEncoding::Quadrant:
            return {2, 2};
        case TextEncoding::Sextant:
            return {2, 3};
        case TextEncoding::Braille:
            return {2, 4};
        default:
            return {1, 1};
    }
}

void tx_utf8(std::string& rOut, std::uint32_t Cp) {
    if (Cp < 0x80) {
        rOut.push_back(static_cast<char>(Cp));
    } else if (Cp < 0x800) {
        rOut.push_back(static_cast<char>(0xC0 | (Cp >> 6)));
        rOut.push_back(static_cast<char>(0x80 | (Cp & 0x3F)));
    } else if (Cp < 0x10000) {
        rOut.push_back(static_cast<char>(0xE0 | (Cp >> 12)));
        rOut.push_back(static_cast<char>(0x80 | ((Cp >> 6) & 0x3F)));
        rOut.push_back(static_cast<char>(0x80 | (Cp & 0x3F)));
    } else {
        rOut.push_back(static_cast<char>(0xF0 | (Cp >> 18)));
        rOut.push_back(static_cast<char>(0x80 | ((Cp >> 12) & 0x3F)));
        rOut.push_back(static_cast<char>(0x80 | ((Cp >> 6) & 0x3F)));
        rOut.push_back(static_cast<char>(0x80 | (Cp & 0x3F)));
    }
}

// The code point that shows the subpixels in Mask (bit i set: subpixel i is
// drawn in the foreground colour).
std::uint32_t tx_glyph(TextEncoding Encoding, std::uint32_t Mask) {
    switch (Encoding) {
        case TextEncoding::HalfBlock: {
            static const std::uint32_t table[4] = {0x20, 0x2580, 0x2584, 0x2588};
            return table[Mask & 3U];
        }
        case TextEncoding::Quadrant: {
            // Bits: 1 upper left, 2 upper right, 4 lower left, 8 lower right.
            static const std::uint32_t table[16] = {0x20,   0x2598, 0x259D, 0x2580, 0x2596, 0x258C,
                                                    0x259E, 0x259B, 0x2597, 0x259A, 0x2590, 0x259C,
                                                    0x2584, 0x2599, 0x259F, 0x2588};
            return table[Mask & 15U];
        }
        case TextEncoding::Sextant: {
            // U+1FB00..U+1FB3B are the sextants in mask order, minus the four that
            // exist elsewhere: empty, left half (1+4+16), right half (2+8+32), full.
            if (Mask == 0)
                return 0x20;
            if (Mask == 63)
                return 0x2588;
            if (Mask == 21)
                return 0x258C;
            if (Mask == 42)
                return 0x2590;
            return 0x1FB00 + Mask - 1 - (Mask > 21 ? 1 : 0) - (Mask > 42 ? 1 : 0);
        }
        case TextEncoding::Braille: {
            // Row-major subpixels to Braille dots 1,4 / 2,5 / 3,6 / 7,8.
            static const std::uint32_t dot[8] = {0x01, 0x08, 0x02, 0x10, 0x04, 0x20, 0x40, 0x80};
            std::uint32_t bits = 0;
            for (int i = 0; i < 8; ++i)
                if (Mask & (1U << i))
                    bits |= dot[i];
            return 0x2800 + bits;
        }
        default:
            return 0x20;
    }
}

// Luma in [0, 255]: BT.709 weights as integers that sum to 256.
int tx_luma(int R, int G, int B) {
    return (54 * R + 183 * G + 19 * B) >> 8;
}

constexpr const char kTxRamp[] = " .:-=+*#%@";

// ---------------------------------------------------------------------------
// Colour depth
// ---------------------------------------------------------------------------

constexpr int kTxCube[6] = {0, 95, 135, 175, 215, 255};

// The xterm defaults for the 16 ANSI colours.
constexpr int kTxAnsi16[16][3] = {{0, 0, 0},       {205, 0, 0},   {0, 205, 0},   {205, 205, 0},
                                  {0, 0, 238},     {205, 0, 205}, {0, 205, 205}, {229, 229, 229},
                                  {127, 127, 127}, {255, 0, 0},   {0, 255, 0},   {255, 255, 0},
                                  {92, 92, 255},   {255, 0, 255}, {0, 255, 255}, {255, 255, 255}};

int tx_dist2(const TxRgb& rA, int R, int G, int B) {
    const int dr = rA[0] - R;
    const int dg = rA[1] - G;
    const int db = rA[2] - B;
    return dr * dr + dg * dg + db * db;
}

int tx_cube_level(int c) {
    int best = 0;
    for (int i = 1; i < 6; ++i)
        if (std::abs(kTxCube[i] - c) < std::abs(kTxCube[best] - c))
            best = i;
    return best;
}

// The nearest xterm-256 entry by squared RGB distance among the nearest cube
// colour and the nearest grey; the cube wins a tie.
int tx_index256(const TxRgb& rC) {
    const int r = tx_cube_level(rC[0]);
    const int g = tx_cube_level(rC[1]);
    const int b = tx_cube_level(rC[2]);
    const int cube = 16 + 36 * r + 6 * g + b;
    const int cube_d = tx_dist2(rC, kTxCube[r], kTxCube[g], kTxCube[b]);
    const int mean = (rC[0] + rC[1] + rC[2]) / 3;
    int grey = (mean - 8 + 5) / 10;
    grey = std::min(std::max(grey, 0), 23);
    const int gv = 8 + 10 * grey;
    const int grey_d = tx_dist2(rC, gv, gv, gv);
    return grey_d < cube_d ? 232 + grey : cube;
}

TxRgb tx_rgb256(int Index) {
    if (Index < 16)
        return {kTxAnsi16[Index][0], kTxAnsi16[Index][1], kTxAnsi16[Index][2]};
    if (Index >= 232) {
        const int v = 8 + 10 * (Index - 232);
        return {v, v, v};
    }
    const int i = Index - 16;
    return {kTxCube[i / 36], kTxCube[(i / 6) % 6], kTxCube[i % 6]};
}

int tx_index16(const TxRgb& rC) {
    int best = 0;
    int best_d = tx_dist2(rC, kTxAnsi16[0][0], kTxAnsi16[0][1], kTxAnsi16[0][2]);
    for (int i = 1; i < 16; ++i) {
        const int d = tx_dist2(rC, kTxAnsi16[i][0], kTxAnsi16[i][1], kTxAnsi16[i][2]);
        if (d < best_d) {
            best = i;
            best_d = d;
        }
    }
    return best;
}

// A colour as the chosen depth shows it.
TxRgb tx_quantize(const TxRgb& rC, ColorDepth Depth) {
    if (Depth == ColorDepth::Palette256)
        return tx_rgb256(tx_index256(rC));
    if (Depth == ColorDepth::Ansi16)
        return tx_rgb256(tx_index16(rC));
    return rC;
}

// The SGR parameters selecting a colour; Fg false for the background. An unset
// colour (the terminal's own) is 39 / 49.
std::string tx_sgr(bool Set, const TxRgb& rC, ColorDepth Depth, bool Fg) {
    if (!Set)
        return Fg ? "39" : "49";
    if (Depth == ColorDepth::Palette256)
        return std::string(Fg ? "38;5;" : "48;5;") + std::to_string(tx_index256(rC));
    if (Depth == ColorDepth::Ansi16) {
        const int i = tx_index16(rC);
        const int base = Fg ? (i < 8 ? 30 : 90 - 8) : (i < 8 ? 40 : 100 - 8);
        return std::to_string(base + i);
    }
    return std::string(Fg ? "38;2;" : "48;2;") + std::to_string(rC[0]) + ";" +
           std::to_string(rC[1]) + ";" + std::to_string(rC[2]);
}

std::string tx_hex(const TxRgb& rC) {
    static const char digits[] = "0123456789abcdef";
    std::string out = "#";
    for (int c : rC) {
        out.push_back(digits[(c >> 4) & 15]);
        out.push_back(digits[c & 15]);
    }
    return out;
}

// ---------------------------------------------------------------------------
// The two-colour split of one cell
// ---------------------------------------------------------------------------

using TxCell = TextCell;  // the public cell: the loop diffs and repaints these

struct TxSub {
    bool mOpaque = false;
    int mC[3] = {0, 0, 0};
};

// Mono: a subpixel is drawn when it is opaque and bright, so dark edges and
// shadows read as gaps in the glyph.
TxCell tx_split_mono(const TxSub* pSub, int Count, TextEncoding Encoding) {
    std::uint32_t mask = 0;
    for (int i = 0; i < Count; ++i)
        if (pSub[i].mOpaque && tx_luma(pSub[i].mC[0], pSub[i].mC[1], pSub[i].mC[2]) >= 128)
            mask |= 1U << i;
    TxCell cell;
    cell.mGlyph = tx_glyph(Encoding, mask);
    return cell;
}

// The split of the subpixels into a foreground and a background group with the
// least total squared error. Transparent subpixels form their own group (cost
// 0) and never share one with opaque ones. Candidates are the masks that put
// subpixel 0 in the foreground: the single group (every subpixel) first, then
// the two-group masks in ascending order; the first minimum wins, so a cell
// that one colour describes exactly never becomes a glyph.
TxCell tx_split_color(const TxSub* pSub, int Count, TextEncoding Encoding) {
    const std::uint32_t full = (1U << Count) - 1U;
    std::uint32_t best_mask = full;
    double best_cost = -1.0;
    std::int64_t best_sum[2][3] = {{0, 0, 0}, {0, 0, 0}};
    std::int64_t best_n[2] = {0, 0};
    bool best_opaque[2] = {false, false};
    for (std::uint32_t step = 0; step <= full / 2; ++step) {
        const std::uint32_t mask = step == 0 ? full : 2 * step - 1;
        std::int64_t n[2] = {0, 0};
        std::int64_t sum[2][3] = {{0, 0, 0}, {0, 0, 0}};
        std::int64_t sum2[2] = {0, 0};
        int opaque[2] = {0, 0};
        for (int i = 0; i < Count; ++i) {
            const int g = (mask >> i) & 1U ? 0 : 1;
            ++n[g];
            if (pSub[i].mOpaque) {
                ++opaque[g];
                for (int c = 0; c < 3; ++c) {
                    sum[g][c] += pSub[i].mC[c];
                    sum2[g] += static_cast<std::int64_t>(pSub[i].mC[c]) * pSub[i].mC[c];
                }
            }
        }
        bool valid = true;
        double cost = 0.0;
        for (int g = 0; g < 2; ++g) {
            if (n[g] == 0)
                continue;
            if (opaque[g] != 0 && opaque[g] != n[g]) {
                valid = false;
                break;
            }
            if (opaque[g] == 0)
                continue;
            const std::int64_t num =
                n[g] * sum2[g] -
                (sum[g][0] * sum[g][0] + sum[g][1] * sum[g][1] + sum[g][2] * sum[g][2]);
            cost += static_cast<double>(num) / static_cast<double>(n[g]);
        }
        if (!valid)
            continue;
        if (best_cost < 0.0 || cost < best_cost) {
            best_cost = cost;
            best_mask = mask;
            for (int g = 0; g < 2; ++g) {
                best_n[g] = n[g];
                best_opaque[g] = opaque[g] > 0;
                for (int c = 0; c < 3; ++c)
                    best_sum[g][c] = sum[g][c];
            }
        }
    }
    TxCell cell;
    TxRgb mean[2];
    for (int g = 0; g < 2; ++g)
        for (int c = 0; c < 3; ++c)
            mean[g][c] =
                best_n[g] > 0 ? static_cast<int>((best_sum[g][c] + best_n[g] / 2) / best_n[g]) : 0;
    std::uint32_t mask = best_mask;
    bool fg_set = best_opaque[0];
    bool bg_set = best_n[1] > 0 && best_opaque[1];
    TxRgb fg = mean[0];
    TxRgb bg = mean[1];
    if (mask == full) {
        // One group: a space on that background (no glyph edge to show).
        cell.mGlyph = 0x20;
        cell.mBgSet = fg_set;
        cell.mBg = fg;
        return cell;
    }
    if (!fg_set && bg_set) {
        // The default foreground is the text colour, not the background: draw
        // the opaque group as the glyph instead.
        mask = full & ~mask;
        std::swap(fg, bg);
        fg_set = true;
        bg_set = false;
    }
    cell.mGlyph = tx_glyph(Encoding, mask);
    cell.mFgSet = fg_set;
    cell.mBgSet = bg_set;
    cell.mFg = fg;
    cell.mBg = bg;
    if (!fg_set && !bg_set)
        cell.mGlyph = 0x20;
    return cell;
}

// The cells of a frame, row by row.
std::vector<TxCell> tx_cells(const Frame& rFrame, const TextOptions& rOpt, int& rCols, int& rRows) {
    const std::array<int, 2> grid = tx_grid(rOpt.mEncoding);
    if (rFrame.mWidth % grid[0] != 0 || rFrame.mHeight % grid[1] != 0)
        throw std::invalid_argument(std::string(kTxPrefix) + "a " + std::to_string(rFrame.mWidth) +
                                    "x" + std::to_string(rFrame.mHeight) +
                                    " frame is not a whole number of cells of this encoding (" +
                                    std::to_string(grid[0]) + "x" + std::to_string(grid[1]) +
                                    " pixels); size it with text_frame_size");
    rCols = rFrame.mWidth / grid[0];
    rRows = rFrame.mHeight / grid[1];
    std::vector<TxCell> cells(static_cast<std::size_t>(rCols) * static_cast<std::size_t>(rRows));
    const bool mono = rOpt.mDepth == ColorDepth::Mono || rOpt.mFormat == TextFormat::Plain;
    for (int row = 0; row < rRows; ++row) {
        for (int col = 0; col < rCols; ++col) {
            TxSub sub[8];
            int count = 0;
            for (int sy = 0; sy < grid[1]; ++sy) {
                for (int sx = 0; sx < grid[0]; ++sx) {
                    const std::size_t px = static_cast<std::size_t>(col * grid[0] + sx);
                    const std::size_t py = static_cast<std::size_t>(row * grid[1] + sy);
                    const std::uint8_t* p = rFrame.mRgba.data() +
                                            (py * static_cast<std::size_t>(rFrame.mWidth) + px) * 4;
                    TxSub& s = sub[count++];
                    s.mOpaque = p[3] >= 128;
                    s.mC[0] = p[0];
                    s.mC[1] = p[1];
                    s.mC[2] = p[2];
                }
            }
            TxCell& cell = cells[static_cast<std::size_t>(row) * static_cast<std::size_t>(rCols) +
                                 static_cast<std::size_t>(col)];
            if (rOpt.mEncoding == TextEncoding::Ascii) {
                if (sub[0].mOpaque) {
                    cell.mGlyph = static_cast<unsigned char>(
                        kTxRamp[tx_luma(sub[0].mC[0], sub[0].mC[1], sub[0].mC[2]) * 10 / 256]);
                    if (!mono) {
                        cell.mFgSet = true;
                        cell.mFg = {sub[0].mC[0], sub[0].mC[1], sub[0].mC[2]};
                    }
                }
            } else if (mono) {
                cell = tx_split_mono(sub, count, rOpt.mEncoding);
            } else {
                cell = tx_split_color(sub, count, rOpt.mEncoding);
            }
        }
    }
    return cells;
}

void tx_html_escape(std::string& rOut, const std::string& rText) {
    for (char c : rText) {
        if (c == '<')
            rOut += "&lt;";
        else if (c == '>')
            rOut += "&gt;";
        else if (c == '&')
            rOut += "&amp;";
        else
            rOut.push_back(c);
    }
}

std::string tx_encode_cells(const Frame& rFrame, const TextOptions& rOpt) {
    int cols = 0;
    int rows = 0;
    const std::vector<TxCell> cells = tx_cells(rFrame, rOpt, cols, rows);
    std::string out;
    const bool html = rOpt.mFormat == TextFormat::Html;
    const bool sgr = rOpt.mFormat == TextFormat::Ansi && rOpt.mDepth != ColorDepth::Mono;
    if (html)
        out +=
            "<!DOCTYPE html>\n<html><head><meta charset=\"utf-8\"><title>meshio++ render</title>\n"
            "<style>pre{font-family:\"DejaVu Sans Mono\",Menlo,Consolas,\"Liberation "
            "Mono\",monospace;"
            "font-size:14px;line-height:1;letter-spacing:0;margin:0}</style>\n</head><body><pre>";
    for (int row = 0; row < rows; ++row) {
        std::string cur_fg = "39";
        std::string cur_bg = "49";
        std::string open_style;
        bool span_open = false;
        for (int col = 0; col < cols; ++col) {
            const TxCell& cell =
                cells[static_cast<std::size_t>(row) * static_cast<std::size_t>(cols) +
                      static_cast<std::size_t>(col)];
            if (sgr) {
                const std::string fg = tx_sgr(cell.mFgSet, cell.mFg, rOpt.mDepth, true);
                const std::string bg = tx_sgr(cell.mBgSet, cell.mBg, rOpt.mDepth, false);
                std::string params;
                if (fg != cur_fg)
                    params = fg;
                if (bg != cur_bg)
                    params += (params.empty() ? "" : ";") + bg;
                if (!params.empty())
                    out += "\x1b[" + params + "m";
                cur_fg = fg;
                cur_bg = bg;
            } else if (html) {
                std::string style;
                if (cell.mFgSet && rOpt.mDepth != ColorDepth::Mono)
                    style += "color:" + tx_hex(tx_quantize(cell.mFg, rOpt.mDepth));
                if (cell.mBgSet && rOpt.mDepth != ColorDepth::Mono)
                    style += std::string(style.empty() ? "" : ";") +
                             "background:" + tx_hex(tx_quantize(cell.mBg, rOpt.mDepth));
                if (!span_open || style != open_style) {
                    if (span_open)
                        out += "</span>";
                    span_open = !style.empty();
                    open_style = style;
                    if (span_open)
                        out += "<span style=\"" + style + "\">";
                }
            }
            tx_utf8(out, cell.mGlyph);
        }
        if (sgr && (cur_fg != "39" || cur_bg != "49"))
            out += "\x1b[0m";
        if (html && span_open)
            out += "</span>";
        out += "\n";
    }
    if (html) {
        out += "</pre>";
        if (rOpt.mNotes)
            for (const std::string& note : rFrame.mNotes) {
                out += "\n<p>";
                tx_html_escape(out, note);
                out += "</p>";
            }
        out += "\n</body></html>\n";
    } else if (rOpt.mNotes) {
        for (const std::string& note : rFrame.mNotes)
            out += note + "\n";
    }
    return out;
}

// ---------------------------------------------------------------------------
// Graphics protocols
// ---------------------------------------------------------------------------

// tmux's passthrough: the sequence inside a DCS, every ESC doubled.
std::string tx_tmux(const std::string& rSeq) {
    std::string out = "\x1bPtmux;";
    for (char c : rSeq) {
        if (c == '\x1b')
            out.push_back('\x1b');
        out.push_back(c);
    }
    out += "\x1b\\";
    return out;
}

std::string tx_kitty(const Frame& rFrame, bool Tmux) {
    const std::string payload = detail::b64encode(rFrame.mRgba.data(), rFrame.mRgba.size());
    std::string out;
    const std::size_t chunk = 4096;
    std::size_t pos = 0;
    bool first = true;
    do {
        const std::size_t len = std::min(chunk, payload.size() - pos);
        const bool more = pos + len < payload.size();
        std::string seq = "\x1b_G";
        if (first)
            seq += "a=T,f=32,s=" + std::to_string(rFrame.mWidth) +
                   ",v=" + std::to_string(rFrame.mHeight) + ",q=2,";
        seq += std::string("m=") + (more ? "1" : "0") + ";" + payload.substr(pos, len) + "\x1b\\";
        out += Tmux ? tx_tmux(seq) : seq;
        pos += len;
        first = false;
    } while (pos < payload.size());
    out += "\n";
    return out;
}

std::string tx_iterm2(const Frame& rFrame, bool Tmux) {
    const std::string png =
        detail::png_encode_rgba(rFrame.mRgba.data(), static_cast<std::size_t>(rFrame.mWidth),
                                static_cast<std::size_t>(rFrame.mHeight), 0);
    const std::string seq =
        "\x1b]1337;File=inline=1;size=" + std::to_string(png.size()) +
        ";width=" + std::to_string(rFrame.mWidth) + "px;height=" + std::to_string(rFrame.mHeight) +
        "px;preserveAspectRatio=1:" +
        detail::b64encode(reinterpret_cast<const unsigned char*>(png.data()), png.size()) + "\a";
    return (Tmux ? tx_tmux(seq) : seq) + "\n";
}

// A deterministic median cut over the opaque colours, at most 256 entries.
struct TxPalette {
    std::vector<std::uint32_t> mColors;  // sorted unique 0xRRGGBB
    std::vector<int> mIndex;             // palette entry per unique colour
    std::vector<TxRgb> mEntries;
};

int tx_channel(std::uint32_t c, int Ch) {
    return static_cast<int>((c >> (16 - 8 * Ch)) & 0xFFU);
}

TxPalette tx_median_cut(const Frame& rFrame) {
    TxPalette pal;
    const std::size_t n =
        static_cast<std::size_t>(rFrame.mWidth) * static_cast<std::size_t>(rFrame.mHeight);
    std::vector<std::uint32_t> all;
    all.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint8_t* p = rFrame.mRgba.data() + i * 4;
        if (p[3] >= 128)
            all.push_back((static_cast<std::uint32_t>(p[0]) << 16) |
                          (static_cast<std::uint32_t>(p[1]) << 8) | p[2]);
    }
    std::sort(all.begin(), all.end());
    std::vector<std::uint64_t> counts;
    for (std::size_t i = 0; i < all.size(); ++i) {
        if (pal.mColors.empty() || pal.mColors.back() != all[i]) {
            pal.mColors.push_back(all[i]);
            counts.push_back(0);
        }
        ++counts.back();
    }
    pal.mIndex.assign(pal.mColors.size(), 0);
    if (pal.mColors.empty())
        return pal;
    // Each box is a list of unique-colour indices.
    std::vector<std::vector<std::size_t>> boxes(1);
    for (std::size_t i = 0; i < pal.mColors.size(); ++i)
        boxes[0].push_back(i);
    while (boxes.size() < 256) {
        int best_box = -1;
        int best_ch = 0;
        int best_range = 0;
        for (std::size_t b = 0; b < boxes.size(); ++b) {
            if (boxes[b].size() < 2)
                continue;
            for (int ch = 0; ch < 3; ++ch) {
                int lo = 255, hi = 0;
                for (std::size_t i : boxes[b]) {
                    lo = std::min(lo, tx_channel(pal.mColors[i], ch));
                    hi = std::max(hi, tx_channel(pal.mColors[i], ch));
                }
                if (hi - lo > best_range) {
                    best_range = hi - lo;
                    best_box = static_cast<int>(b);
                    best_ch = ch;
                }
            }
        }
        if (best_box < 0)
            break;
        std::vector<std::size_t>& box = boxes[static_cast<std::size_t>(best_box)];
        std::sort(box.begin(), box.end(), [&](std::size_t a, std::size_t b) {
            const int ca = tx_channel(pal.mColors[a], best_ch);
            const int cb = tx_channel(pal.mColors[b], best_ch);
            return ca != cb ? ca < cb : pal.mColors[a] < pal.mColors[b];
        });
        std::uint64_t total = 0;
        for (std::size_t i : box)
            total += counts[i];
        std::uint64_t acc = 0;
        std::size_t cut = 1;
        for (std::size_t k = 0; k + 1 < box.size(); ++k) {
            acc += counts[box[k]];
            cut = k + 1;
            if (2 * acc >= total)
                break;
        }
        std::vector<std::size_t> upper(box.begin() + static_cast<std::ptrdiff_t>(cut), box.end());
        box.resize(cut);
        boxes.push_back(std::move(upper));
    }
    for (std::size_t b = 0; b < boxes.size(); ++b) {
        std::uint64_t w = 0;
        std::uint64_t sum[3] = {0, 0, 0};
        for (std::size_t i : boxes[b]) {
            w += counts[i];
            for (int ch = 0; ch < 3; ++ch)
                sum[ch] += counts[i] * static_cast<std::uint64_t>(tx_channel(pal.mColors[i], ch));
            pal.mIndex[i] = static_cast<int>(b);
        }
        TxRgb e;
        for (int ch = 0; ch < 3; ++ch)
            e[ch] = static_cast<int>((sum[ch] + w / 2) / w);
        pal.mEntries.push_back(e);
    }
    return pal;
}

std::string tx_sixel(const Frame& rFrame, bool Tmux) {
    const TxPalette pal = tx_median_cut(rFrame);
    const int w = rFrame.mWidth;
    const int h = rFrame.mHeight;
    // The palette index of each pixel, -1 for a transparent one.
    std::vector<int> px(static_cast<std::size_t>(w) * static_cast<std::size_t>(h), -1);
    for (std::size_t i = 0; i < px.size(); ++i) {
        const std::uint8_t* p = rFrame.mRgba.data() + i * 4;
        if (p[3] < 128)
            continue;
        const std::uint32_t c = (static_cast<std::uint32_t>(p[0]) << 16) |
                                (static_cast<std::uint32_t>(p[1]) << 8) | p[2];
        const auto it = std::lower_bound(pal.mColors.begin(), pal.mColors.end(), c);
        px[i] = pal.mIndex[static_cast<std::size_t>(it - pal.mColors.begin())];
    }
    std::string seq = "\x1bP0;1;0q\"1;1;" + std::to_string(w) + ";" + std::to_string(h);
    for (std::size_t e = 0; e < pal.mEntries.size(); ++e) {
        seq += "#" + std::to_string(e) + ";2";
        for (int ch = 0; ch < 3; ++ch)
            seq += ";" + std::to_string((pal.mEntries[e][ch] * 100 + 127) / 255);
    }
    auto run = [&](std::string& rOut, char Ch, int Count) {
        if (Count >= 4)
            rOut += "!" + std::to_string(Count) + Ch;
        else
            rOut.append(static_cast<std::size_t>(Count), Ch);
    };
    for (int band = 0; band < h; band += 6) {
        bool first_color = true;
        for (std::size_t e = 0; e < pal.mEntries.size(); ++e) {
            std::string line;
            bool used = false;
            char prev = 0;
            int count = 0;
            for (int x = 0; x < w; ++x) {
                int bits = 0;
                for (int r = 0; r < 6 && band + r < h; ++r)
                    if (px[static_cast<std::size_t>(band + r) * static_cast<std::size_t>(w) +
                           static_cast<std::size_t>(x)] == static_cast<int>(e))
                        bits |= 1 << r;
                if (bits)
                    used = true;
                const char ch = static_cast<char>(63 + bits);
                if (count > 0 && ch == prev) {
                    ++count;
                } else {
                    if (count > 0)
                        run(line, prev, count);
                    prev = ch;
                    count = 1;
                }
            }
            if (!used)
                continue;
            if (count > 0 && prev != 63)
                run(line, prev, count);
            seq += (first_color ? "" : "$") + std::string("#") + std::to_string(e) + line;
            first_color = false;
        }
        seq += "-";
    }
    seq += "\x1b\\";
    return (Tmux ? tx_tmux(seq) : seq) + "\n";
}

// ---------------------------------------------------------------------------
// asciicast
// ---------------------------------------------------------------------------

void tx_json_string(std::string& rOut, const std::string& rText) {
    rOut.push_back('"');
    for (unsigned char c : rText) {
        if (c == '"') {
            rOut += "\\\"";
        } else if (c == '\\') {
            rOut += "\\\\";
        } else if (c < 0x20) {
            char buf[8];
            detail::snprintf_c(buf, sizeof(buf), "\\u%04x", static_cast<unsigned>(c));
            rOut += buf;
        } else {
            rOut.push_back(static_cast<char>(c));
        }
    }
    rOut.push_back('"');
}

std::string tx_lower_ext(const std::string& rPath) {
    const std::size_t slash = rPath.find_last_of("/\\");
    const std::size_t dot = rPath.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
        return "";
    std::string ext = rPath.substr(dot);
    for (char& c : ext)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext;
}

}  // namespace

std::array<int, 2> text_frame_size(const TextOptions& rOptions, double& rPixelAspect) {
    if (rOptions.mCols <= 0 || rOptions.mRows <= 0)
        throw std::invalid_argument(std::string(kTxPrefix) + "columns and rows must be positive");
    if (!(rOptions.mCellAspect > 0.0) || !std::isfinite(rOptions.mCellAspect))
        throw std::invalid_argument(std::string(kTxPrefix) + "cell aspect must be positive");
    if (!tx_is_cell(rOptions.mEncoding)) {
        if (rOptions.mCellPixelWidth <= 0 || rOptions.mCellPixelHeight <= 0)
            throw std::invalid_argument(std::string(kTxPrefix) +
                                        "cell pixel size must be positive");
        rPixelAspect = 1.0;
        return {rOptions.mCols * rOptions.mCellPixelWidth,
                rOptions.mRows * rOptions.mCellPixelHeight};
    }
    const std::array<int, 2> grid = tx_grid(rOptions.mEncoding);
    rPixelAspect =
        rOptions.mCellAspect * static_cast<double>(grid[0]) / static_cast<double>(grid[1]);
    return {rOptions.mCols * grid[0], rOptions.mRows * grid[1]};
}

std::string encode_text(const Frame& rFrame, const TextOptions& rOptions) {
    if (rFrame.mWidth <= 0 || rFrame.mHeight <= 0 ||
        rFrame.mRgba.size() !=
            static_cast<std::size_t>(rFrame.mWidth) * static_cast<std::size_t>(rFrame.mHeight) * 4)
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "the frame's buffer does not match its size");
    if (tx_is_cell(rOptions.mEncoding))
        return tx_encode_cells(rFrame, rOptions);
    if (rOptions.mFormat != TextFormat::Ansi)
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "the Kitty, iTerm2 and Sixel protocols are escape sequences; "
                                    "they need the ANSI format, not plain text or HTML");
    switch (rOptions.mEncoding) {
        case TextEncoding::Kitty:
            return tx_kitty(rFrame, rOptions.mTmuxPassthrough);
        case TextEncoding::ITerm2:
            return tx_iterm2(rFrame, rOptions.mTmuxPassthrough);
        default:
            return tx_sixel(rFrame, rOptions.mTmuxPassthrough);
    }
}

TextGrid encode_cells(const Frame& rFrame, const TextOptions& rOptions) {
    if (rFrame.mWidth <= 0 || rFrame.mHeight <= 0 ||
        rFrame.mRgba.size() !=
            static_cast<std::size_t>(rFrame.mWidth) * static_cast<std::size_t>(rFrame.mHeight) * 4)
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "the frame's buffer does not match its size");
    if (!tx_is_cell(rOptions.mEncoding))
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "encode_cells needs a cell encoding, not a graphics protocol");
    TextGrid grid;
    grid.mCells = tx_cells(rFrame, rOptions, grid.mCols, grid.mRows);
    return grid;
}

namespace {

// Whether two cells look the same at a colour depth.
bool tx_same_cell(const TextCell& rA, const TextCell& rB, ColorDepth Depth) {
    if (rA.mGlyph != rB.mGlyph)
        return false;
    if (Depth == ColorDepth::Mono)
        return true;
    if (rA.mFgSet != rB.mFgSet || rA.mBgSet != rB.mBgSet)
        return false;
    if (rA.mFgSet && tx_quantize(rA.mFg, Depth) != tx_quantize(rB.mFg, Depth))
        return false;
    if (rA.mBgSet && tx_quantize(rA.mBg, Depth) != tx_quantize(rB.mBg, Depth))
        return false;
    return true;
}

}  // namespace

std::string encode_cells_update(const TextGrid& rOld, const TextGrid& rNew, ColorDepth Depth,
                                int Row, int Col) {
    const bool full = rOld.mCols != rNew.mCols || rOld.mRows != rNew.mRows ||
                      rOld.mCells.size() != rNew.mCells.size();
    const bool sgr = Depth != ColorDepth::Mono;
    std::string out;
    std::string cur_fg;  // empty: not known, so the next cell sets both
    std::string cur_bg;
    bool any = false;
    for (int row = 0; row < rNew.mRows; ++row) {
        int last = -2;  // the column the cursor is already at, after the last glyph written
        for (int col = 0; col < rNew.mCols; ++col) {
            const std::size_t at = static_cast<std::size_t>(row) *
                                       static_cast<std::size_t>(rNew.mCols) +
                                   static_cast<std::size_t>(col);
            const TextCell& cell = rNew.mCells[at];
            if (!full && tx_same_cell(rOld.mCells[at], cell, Depth))
                continue;
            if (col != last + 1)
                out += "\x1b[" + std::to_string(Row + row) + ";" + std::to_string(Col + col) + "H";
            if (sgr) {
                const std::string fg = tx_sgr(cell.mFgSet, cell.mFg, Depth, true);
                const std::string bg = tx_sgr(cell.mBgSet, cell.mBg, Depth, false);
                std::string params;
                if (fg != cur_fg)
                    params = fg;
                if (bg != cur_bg)
                    params += (params.empty() ? "" : ";") + bg;
                if (!params.empty())
                    out += "\x1b[" + params + "m";
                cur_fg = fg;
                cur_bg = bg;
            }
            tx_utf8(out, cell.mGlyph);
            last = col;
            any = true;
        }
    }
    if (any && sgr)
        out += "\x1b[0m";
    return out;
}

std::string render_text(const Mesh& rMesh, const RenderOptions& rRender, const TextOptions& rText) {
    RenderOptions options = rRender;
    const std::array<int, 2> size = text_frame_size(rText, options.mPixelAspect);
    options.mWidth = size[0];
    options.mHeight = size[1];
    return encode_text(render(rMesh, options), rText);
}

std::string encode_png(const Frame& rFrame, int CompressLevel) {
    if (rFrame.mWidth <= 0 || rFrame.mHeight <= 0 ||
        rFrame.mRgba.size() !=
            static_cast<std::size_t>(rFrame.mWidth) * static_cast<std::size_t>(rFrame.mHeight) * 4)
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "the frame's buffer does not match its size");
    return detail::png_encode_rgba(rFrame.mRgba.data(), static_cast<std::size_t>(rFrame.mWidth),
                                   static_cast<std::size_t>(rFrame.mHeight), CompressLevel);
}

ColorDepth detect_color_depth(const char* pNoColor, const char* pColorTerm, const char* pTerm) {
    if (pNoColor != nullptr && pNoColor[0] != '\0')
        return ColorDepth::Mono;
    if (pColorTerm != nullptr &&
        (std::strcmp(pColorTerm, "truecolor") == 0 || std::strcmp(pColorTerm, "24bit") == 0))
        return ColorDepth::TrueColor;
    if (pTerm != nullptr) {
        if (std::strstr(pTerm, "256color") != nullptr)
            return ColorDepth::Palette256;
        if (std::strcmp(pTerm, "dumb") == 0)
            return ColorDepth::Mono;
    }
    return ColorDepth::Ansi16;
}

std::string encode_cast(const Mesh& rMesh, const RenderOptions& rRender, const TextOptions& rText,
                        const SnapshotOptions& rSnapshot) {
    if (!tx_is_cell(rText.mEncoding))
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "an asciicast records text cells; "
                                    "choose a cell encoding, not a graphics protocol");
    if (rSnapshot.mCastFrames <= 0)
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "an asciicast needs at least one frame");
    if (!(rSnapshot.mCastFps > 0.0) || !std::isfinite(rSnapshot.mCastFps) ||
        !std::isfinite(rSnapshot.mCastDegrees))
        throw std::invalid_argument(std::string(kTxPrefix) +
                                    "an asciicast needs a positive, finite frame rate");
    TextOptions text = rText;
    text.mFormat = TextFormat::Ansi;
    RenderOptions options = rRender;
    // Orbit from the named view, if any: resolve it once, then sweep the azimuth.
    detail::render_view_angles(options, options.mAzimuth, options.mElevation);
    options.mView.clear();
    const double azimuth = options.mAzimuth;
    std::vector<std::string> frames;
    int height = text.mRows;
    for (int i = 0; i < rSnapshot.mCastFrames; ++i) {
        options.mAzimuth = azimuth + rSnapshot.mCastDegrees * static_cast<double>(i) /
                                         static_cast<double>(rSnapshot.mCastFrames);
        std::string body = render_text(rMesh, options, text);
        std::string crlf;
        int lines = 0;
        for (char c : body) {
            if (c == '\n') {
                crlf += "\r\n";
                ++lines;
            } else {
                crlf.push_back(c);
            }
        }
        height = std::max(height, lines);
        frames.push_back(std::move(crlf));
    }
    std::string out = "{\"height\": " + std::to_string(height) +
                      ", \"version\": 2, \"width\": " + std::to_string(text.mCols) + "}\n";
    for (std::size_t i = 0; i < frames.size(); ++i) {
        char stamp[64];
        detail::snprintf_c(stamp, sizeof(stamp), "%.6f",
                           static_cast<double>(i) / rSnapshot.mCastFps);
        std::string data =
            (i == 0 ? std::string("\x1b[?25l\x1b[2J") : std::string()) + "\x1b[H" + frames[i];
        if (i + 1 == frames.size())
            data += "\x1b[?25h";
        out += "[" + std::string(stamp) + ", \"o\", ";
        tx_json_string(out, data);
        out += "]\n";
    }
    return out;
}

void write_snapshot(const std::string& rPath, const Mesh& rMesh, const RenderOptions& rRender,
                    const TextOptions& rText, const SnapshotOptions& rSnapshot) {
    const std::string ext = tx_lower_ext(rPath);
    std::string bytes;
    TextOptions text = rText;
    if (ext == ".png") {
        bytes = encode_png(render(rMesh, rRender), rSnapshot.mPngCompress);
    } else if (ext == ".txt" || ext == ".ansi" || ext == ".html" || ext == ".htm") {
        text.mFormat = ext == ".txt" ? TextFormat::Plain
                                     : (ext == ".ansi" ? TextFormat::Ansi : TextFormat::Html);
        if (!tx_is_cell(text.mEncoding))
            throw std::invalid_argument(
                std::string(kTxPrefix) + "'" + ext +
                "' holds text cells; choose a cell encoding, not a graphics protocol");
        bytes = render_text(rMesh, rRender, text);
    } else if (ext == ".cast") {
        bytes = encode_cast(rMesh, rRender, text, rSnapshot);
    } else {
        throw std::invalid_argument(std::string(kTxPrefix) + "cannot tell the snapshot form of '" +
                                    rPath + "' (expected .png, .txt, .ansi, .html or .cast)");
    }
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("meshio++: snapshot: cannot open '" + rPath + "' for writing");
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!out)
        throw WriteError("meshio++: snapshot: failed writing '" + rPath + "'");
}

}  // namespace meshioplusplus
