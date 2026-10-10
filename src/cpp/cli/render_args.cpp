// SPDX-License-Identifier: MIT
#include "render_args.hpp"

#include <cstdint>
#include <stdexcept>

#include "meshioplusplus/detail/fast_number.hpp"
#include "terminal.hpp"

namespace meshioplusplus::cli {

/// Parse `args` against `specs`. Throws std::runtime_error on an unknown option
/// or a value-option missing its argument. `--` stops option parsing.
cli_parsed cli_parse(const std::vector<std::string>& rArgs,
                     const std::vector<cli_opt_spec>& rSpecs) {
    cli_parsed out;
    bool positional_only = false;
    for (std::size_t i = 0; i < rArgs.size(); ++i) {
        const std::string& tok = rArgs[i];
        if (positional_only || tok.empty() || tok[0] != '-' || tok == "-") {
            out.positionals.push_back(tok);
            continue;
        }
        if (tok == "--") {
            positional_only = true;
            continue;
        }
        // Split "--name=value".
        std::string name = tok;
        std::string inline_value;
        bool has_inline = false;
        auto eq = tok.find('=');
        if (tok.rfind("--", 0) == 0 && eq != std::string::npos) {
            name = tok.substr(0, eq);
            inline_value = tok.substr(eq + 1);
            has_inline = true;
        }
        const cli_opt_spec* match = nullptr;
        for (const auto& s : rSpecs) {
            if (name == "--" + s.canonical) {
                match = &s;
                break;
            }
            for (const auto& a : s.aliases) {
                if (name == a) {
                    match = &s;
                    break;
                }
            }
            if (match)
                break;
        }
        if (!match)
            throw std::runtime_error("unknown option '" + name + "'");
        if (match->takes_value) {
            std::string value;
            if (has_inline) {
                value = inline_value;
            } else if (i + 1 < rArgs.size()) {
                value = rArgs[++i];
            } else {
                throw std::runtime_error("option '" + name + "' requires a value");
            }
            out.values[match->canonical] = value;
            out.multi[match->canonical].push_back(value);
            out.ordered.emplace_back(match->canonical, value);
        } else {
            out.flags.insert(match->canonical);
        }
    }
    return out;
}

std::string opt_value(const cli_parsed& rP, const std::string& rName,
                      const std::string& rDefault) {
    auto it = rP.values.find(rName);
    return it == rP.values.end() ? rDefault : it->second;
}
bool has_flag(const cli_parsed& rP, const std::string& rName) {
    return rP.flags.count(rName) != 0;
}
/// Every occurrence of a repeatable value-option, in order (empty if absent).
const std::vector<std::string>& opt_values(const cli_parsed& rP, const std::string& rName) {
    static const std::vector<std::string> empty;
    auto it = rP.multi.find(rName);
    return it == rP.multi.end() ? empty : it->second;
}
/// Whether a value-option was supplied at all (distinct from "supplied empty").
bool has_opt(const cli_parsed& rP, const std::string& rName) {
    return rP.values.count(rName) != 0;
}

/// The numeric half of the `--color-by` family (`--component`, `--vmin`,
/// `--vmax`), shared by `convert` and `snapshot` so both read them alike.
void cli_color_values(const cli_parsed& rP, std::optional<int>& rComponent,
                      std::optional<double>& rVMin, std::optional<double>& rVMax) {
    if (has_opt(rP, "component"))
        rComponent = std::stoi(opt_value(rP, "component"));
    if (has_opt(rP, "vmin"))
        rVMin = meshioplusplus::detail::stod_c(opt_value(rP, "vmin"));
    if (has_opt(rP, "vmax"))
        rVMax = meshioplusplus::detail::stod_c(opt_value(rP, "vmax"));
}

/// The flags every software-rendering verb shares (snapshot now, the TUI's
/// `:` commands later): one vocabulary, parsed by one function.
std::vector<cli_opt_spec> render_flag_specs() {
    return {
        {"view", {}, true},           {"azimuth", {}, true},        {"elevation", {}, true},
        {"roll", {}, true},           {"perspective", {}, false},   {"fov", {}, true},
        {"zoom", {}, true},           {"pan-x", {}, true},          {"pan-y", {}, true},
        {"shading", {}, true},        {"one-sided", {}, false},     {"ambient", {}, true},
        {"split-angle", {}, true},    {"edges", {}, true},          {"feature-angle", {}, true},
        {"edge-color", {}, true},     {"fill", {}, true},           {"line-color", {}, true},
        {"background", {}, true},     {"point-radius", {}, true},   {"supersample", {}, true},
        {"axes", {}, false},          {"scale-bar", {}, false},     {"color-by", {}, true},
        {"component", {}, true},      {"cmap", {}, true},           {"vmin", {}, true},
        {"vmax", {}, true},           {"nan-color", {}, true},      {"colorbar", {}, false},
        {"reduce", {}, true},         {"expr", {}, true},           {"clip", {}, true},
        {"symmetric", {}, false},     {"scale", {}, true},          {"scale-threshold", {}, true},
        {"categorical", {}, false},   {"color-regions", {}, false}, {"category-edges", {}, false},
        {"isolines", {}, true},       {"iso-levels", {}, true},     {"iso-color", {}, true},
        {"vectors", {}, true},        {"vector-count", {}, true},   {"vector-length", {}, true},
        {"vector-color", {}, true},   {"warp", {}, true},           {"warp-scale", {}, true},
        {"warp-outline", {}, false},  {"outline-color", {}, true},  {"diagnostic", {}, true},
        {"quality-metric", {}, true},
    };
}

/// A comma-separated list of numbers; an empty entry is `nullopt` (so
/// `--clip 2,` leaves the high end alone).
std::vector<std::optional<double>> cli_parse_numbers(const std::string& rText, const char* pFlag) {
    std::vector<std::optional<double>> out;
    std::size_t at = 0;
    while (true) {
        const std::size_t comma = rText.find(',', at);
        const std::string item = rText.substr(at, comma == std::string::npos ? comma : comma - at);
        if (item.empty()) {
            out.emplace_back(std::nullopt);
        } else {
            const char* end = nullptr;
            const double v = meshioplusplus::detail::parse_double(item.c_str(), end);
            if (end == item.c_str() || end != item.c_str() + item.size())
                throw std::invalid_argument(std::string("--") + pFlag +
                                            " expects numbers separated by commas, not '" + rText +
                                            "'");
            out.emplace_back(v);
        }
        if (comma == std::string::npos)
            break;
        at = comma + 1;
    }
    return out;
}

/// `#rrggbb`, `#rrggbbaa`, or `none`/`transparent` (alpha 0).
meshioplusplus::RenderColor cli_parse_rgba(const std::string& rText, const char* pFlag) {
    if (rText == "none" || rText == "transparent")
        return {0, 0, 0, 0};
    auto hex = [&](std::size_t At) {
        int v = 0;
        for (std::size_t k = At; k < At + 2; ++k) {
            const char c = rText[k];
            v = v * 16 + (c >= '0' && c <= '9'   ? c - '0'
                          : c >= 'a' && c <= 'f' ? c - 'a' + 10
                          : c >= 'A' && c <= 'F' ? c - 'A' + 10
                                                 : -1000);
        }
        if (v < 0)
            throw std::invalid_argument(std::string("--") + pFlag +
                                        " expects #rrggbb, #rrggbbaa or none, not '" + rText + "'");
        return static_cast<std::uint8_t>(v);
    };
    if ((rText.size() != 7 && rText.size() != 9) || rText[0] != '#')
        throw std::invalid_argument(std::string("--") + pFlag +
                                    " expects #rrggbb, #rrggbbaa or none, not '" + rText + "'");
    return {hex(1), hex(3), hex(5), rText.size() == 9 ? hex(7) : static_cast<std::uint8_t>(255)};
}

meshioplusplus::RenderOptions cli_render_options(const cli_parsed& rP) {
    using meshioplusplus::detail::stod_c;
    meshioplusplus::RenderOptions o;
    o.mView = opt_value(rP, "view");
    if (has_opt(rP, "azimuth"))
        o.mAzimuth = stod_c(opt_value(rP, "azimuth"));
    if (has_opt(rP, "elevation"))
        o.mElevation = stod_c(opt_value(rP, "elevation"));
    if (has_opt(rP, "roll"))
        o.mRoll = stod_c(opt_value(rP, "roll"));
    if (!o.mView.empty() && (has_opt(rP, "azimuth") || has_opt(rP, "elevation")))
        throw std::invalid_argument("--view and --azimuth/--elevation are mutually exclusive");
    if (has_flag(rP, "perspective"))
        o.mProjection = meshioplusplus::RenderProjection::Perspective;
    else if (has_opt(rP, "fov"))
        throw std::invalid_argument("--fov requires --perspective");
    if (has_opt(rP, "fov"))
        o.mFovDeg = stod_c(opt_value(rP, "fov"));
    if (has_opt(rP, "zoom"))
        o.mZoom = stod_c(opt_value(rP, "zoom"));
    if (has_opt(rP, "pan-x"))
        o.mPanX = stod_c(opt_value(rP, "pan-x"));
    if (has_opt(rP, "pan-y"))
        o.mPanY = stod_c(opt_value(rP, "pan-y"));
    const std::string shading = opt_value(rP, "shading", "flat");
    if (shading == "none")
        o.mShading = meshioplusplus::RenderShading::None;
    else if (shading == "flat")
        o.mShading = meshioplusplus::RenderShading::Flat;
    else if (shading == "smooth")
        o.mShading = meshioplusplus::RenderShading::Smooth;
    else
        throw std::invalid_argument("--shading expects none, flat or smooth, not '" + shading +
                                    "'");
    o.mTwoSided = !has_flag(rP, "one-sided");
    if (has_opt(rP, "ambient"))
        o.mAmbient = stod_c(opt_value(rP, "ambient"));
    if (has_opt(rP, "split-angle"))
        o.mSplitAngle = stod_c(opt_value(rP, "split-angle"));
    const std::string edges = opt_value(rP, "edges", "none");
    if (edges == "none")
        o.mEdges = meshioplusplus::RenderEdges::None;
    else if (edges == "all")
        o.mEdges = meshioplusplus::RenderEdges::All;
    else if (edges == "feature")
        o.mEdges = meshioplusplus::RenderEdges::Feature;
    else
        throw std::invalid_argument("--edges expects none, all or feature, not '" + edges + "'");
    if (has_opt(rP, "feature-angle"))
        o.mFeatureAngle = stod_c(opt_value(rP, "feature-angle"));
    if (has_opt(rP, "edge-color"))
        o.mEdgeColor = cli_parse_rgba(opt_value(rP, "edge-color"), "edge-color");
    if (has_opt(rP, "fill"))
        o.mFillColor = cli_parse_rgba(opt_value(rP, "fill"), "fill");
    if (has_opt(rP, "line-color"))
        o.mLineColor = cli_parse_rgba(opt_value(rP, "line-color"), "line-color");
    if (has_opt(rP, "background"))
        o.mBackground = cli_parse_rgba(opt_value(rP, "background"), "background");
    if (has_opt(rP, "point-radius"))
        o.mPointRadius = stod_c(opt_value(rP, "point-radius"));
    if (has_opt(rP, "supersample"))
        o.mSupersample = std::stoi(opt_value(rP, "supersample"));
    o.mAxes = has_flag(rP, "axes");
    o.mScaleBar = has_flag(rP, "scale-bar");
    // Field rendering. A mapped field comes from an array, an expression, or
    // a diagnostic that computes one; only then do the range flags apply.
    o.mReduce = opt_value(rP, "reduce");
    o.mExpr = opt_value(rP, "expr");
    if (has_opt(rP, "clip")) {
        const auto clip = cli_parse_numbers(opt_value(rP, "clip"), "clip");
        if (clip.size() != 2)
            throw std::invalid_argument(
                "--clip expects LOW,HIGH percentiles (either may be empty)");
        o.mClipLow = clip[0];
        o.mClipHigh = clip[1];
    }
    o.mSymmetric = has_flag(rP, "symmetric");
    const std::string scale = opt_value(rP, "scale", "linear");
    if (scale == "linear")
        o.mScale = meshioplusplus::RenderScale::Linear;
    else if (scale == "log")
        o.mScale = meshioplusplus::RenderScale::Log;
    else if (scale == "symlog")
        o.mScale = meshioplusplus::RenderScale::Symlog;
    else
        throw std::invalid_argument("--scale expects linear, log or symlog, not '" + scale + "'");
    if (has_opt(rP, "scale-threshold"))
        o.mScaleThreshold = stod_c(opt_value(rP, "scale-threshold"));
    o.mCategorical = has_flag(rP, "categorical");
    o.mColorRegions = has_flag(rP, "color-regions");
    o.mCategoryEdges = has_flag(rP, "category-edges");
    if (has_opt(rP, "isolines"))
        o.mIsolines = std::stoi(opt_value(rP, "isolines"));
    if (has_opt(rP, "iso-levels"))
        for (const auto& level : cli_parse_numbers(opt_value(rP, "iso-levels"), "iso-levels")) {
            if (!level.has_value())
                throw std::invalid_argument("--iso-levels expects numbers separated by commas");
            o.mIsoLevels.push_back(*level);
        }
    if (has_opt(rP, "iso-color"))
        o.mIsoColor = cli_parse_rgba(opt_value(rP, "iso-color"), "iso-color");
    o.mVectors = opt_value(rP, "vectors");
    if (has_opt(rP, "vector-count"))
        o.mVectorCount = std::stoi(opt_value(rP, "vector-count"));
    if (has_opt(rP, "vector-length"))
        o.mVectorLength = stod_c(opt_value(rP, "vector-length"));
    if (has_opt(rP, "vector-color"))
        o.mVectorColor = cli_parse_rgba(opt_value(rP, "vector-color"), "vector-color");
    o.mWarp = opt_value(rP, "warp");
    if (has_opt(rP, "warp-scale"))
        o.mWarpScale = stod_c(opt_value(rP, "warp-scale"));
    o.mWarpOutline = has_flag(rP, "warp-outline");
    if (has_opt(rP, "outline-color"))
        o.mOutlineColor = cli_parse_rgba(opt_value(rP, "outline-color"), "outline-color");
    const std::string diagnostic = opt_value(rP, "diagnostic", "none");
    using meshioplusplus::RenderDiagnostic;
    static const std::pair<const char*, RenderDiagnostic> diagnostics[] = {
        {"none", RenderDiagnostic::None},
        {"quality", RenderDiagnostic::Quality},
        {"inverted", RenderDiagnostic::Inverted},
        {"degenerate", RenderDiagnostic::Degenerate},
        {"orientation", RenderDiagnostic::Orientation},
        {"free-edges", RenderDiagnostic::FreeEdges},
        {"edge-length", RenderDiagnostic::EdgeLength}};
    bool known_diagnostic = false;
    for (const auto& d : diagnostics)
        if (diagnostic == d.first) {
            o.mDiagnostic = d.second;
            known_diagnostic = true;
        }
    if (!known_diagnostic)
        throw std::invalid_argument(
            "--diagnostic expects none, quality, inverted, degenerate, orientation, "
            "free-edges or edge-length, not '" +
            diagnostic + "'");
    o.mQualityMetric = opt_value(rP, "quality-metric");
    if (has_opt(rP, "quality-metric") && o.mDiagnostic != RenderDiagnostic::Quality)
        throw std::runtime_error("--quality-metric requires --diagnostic quality");

    const bool color = has_opt(rP, "color-by");
    const bool mapped = color || !o.mExpr.empty() || o.mDiagnostic == RenderDiagnostic::Quality ||
                        o.mDiagnostic == RenderDiagnostic::EdgeLength;
    if (!mapped) {
        for (const char* flag :
             {"component", "cmap", "vmin", "vmax", "nan-color", "clip", "scale", "scale-threshold"})
            if (has_opt(rP, flag))
                throw std::runtime_error(
                    std::string("--") + flag +
                    " requires --color-by, --expr or a quality or edge-length diagnostic");
        if (has_flag(rP, "colorbar") || has_flag(rP, "symmetric"))
            throw std::runtime_error(
                "--colorbar and --symmetric require --color-by, --expr or a quality or "
                "edge-length diagnostic");
    } else {
        o.mColorBy = opt_value(rP, "color-by");
        cli_color_values(rP, o.mComponent, o.mVMin, o.mVMax);
        o.mCmap = opt_value(rP, "cmap", "viridis");
        if (has_opt(rP, "nan-color"))
            o.mNanColor = cli_parse_rgba(opt_value(rP, "nan-color"), "nan-color");
        o.mColorbar = has_flag(rP, "colorbar");
    }
    return o;
}

meshioplusplus::TextEncoding cli_text_encoding(const std::string& rName) {
    using meshioplusplus::TextEncoding;
    if (rName == "halfblock")
        return TextEncoding::HalfBlock;
    if (rName == "quadrant")
        return TextEncoding::Quadrant;
    if (rName == "sextant")
        return TextEncoding::Sextant;
    if (rName == "braille")
        return TextEncoding::Braille;
    if (rName == "ascii")
        return TextEncoding::Ascii;
    if (rName == "kitty")
        return TextEncoding::Kitty;
    if (rName == "iterm2")
        return TextEncoding::ITerm2;
    if (rName == "sixel")
        return TextEncoding::Sixel;
    throw std::invalid_argument(
        "--encoding expects halfblock, quadrant, sextant, braille, ascii, kitty, "
        "iterm2 or sixel, not '" +
        rName + "'");
}

meshioplusplus::ColorDepth cli_color_depth(const std::string& rName) {
    using meshioplusplus::ColorDepth;
    if (rName == "auto")
        return meshioplusplus::cli::environment_color_depth();
    if (rName == "truecolor" || rName == "24bit")
        return ColorDepth::TrueColor;
    if (rName == "256")
        return ColorDepth::Palette256;
    if (rName == "16")
        return ColorDepth::Ansi16;
    if (rName == "mono")
        return ColorDepth::Mono;
    throw std::invalid_argument("--color-depth expects auto, truecolor, 256, 16 or mono, not '" +
                                rName + "'");
}

}  // namespace meshioplusplus::cli
