// SPDX-License-Identifier: MIT
#include "render_args.hpp"

#include <algorithm>
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
        {"quality-metric", {}, true}, {"cutaway", {}, true},    {"cutaway-tint", {}, true},
        {"streamlines", {}, true},    {"stream-seeds", {}, true},   {"stream-length", {}, true},
        {"stream-color", {}, true},   {"theme", {}, true},          {"scanlines", {}, false},
        {"bloom", {}, false},         {"fringe", {}, false},        {"grid-floor", {}, false},
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

meshioplusplus::RenderCutaway cli_parse_cutaway(const std::string& rText) {
    RenderCutaway cut;
    const std::size_t colon = rText.find(':');
    if (colon != std::string::npos) {
        const std::string axis = rText.substr(0, colon);
        if (axis.size() != 2 || (axis[0] != '+' && axis[0] != '-') ||
            (axis[1] != 'x' && axis[1] != 'y' && axis[1] != 'z'))
            throw std::invalid_argument("--cutaway expects AXIS:OFFSET with AXIS one of +x -x +y -y "
                                        "+z -z, or six numbers PX,PY,PZ,NX,NY,NZ, not '" +
                                        rText + "'");
        const auto offset = cli_parse_numbers(rText.substr(colon + 1), "cutaway");
        if (offset.size() != 1 || !offset[0].has_value())
            throw std::invalid_argument("--cutaway expects one number after the colon, not '" +
                                        rText + "'");
        const std::size_t k = static_cast<std::size_t>(axis[1] - 'x');
        cut.mPoint = {0.0, 0.0, 0.0};
        cut.mPoint[k] = *offset[0];
        cut.mNormal = {0.0, 0.0, 0.0};
        cut.mNormal[k] = axis[0] == '+' ? 1.0 : -1.0;
        return cut;
    }
    const auto numbers = cli_parse_numbers(rText, "cutaway");
    if (numbers.size() != 6)
        throw std::invalid_argument("--cutaway expects AXIS:OFFSET or six numbers "
                                    "PX,PY,PZ,NX,NY,NZ, not '" +
                                    rText + "'");
    for (std::size_t k = 0; k < 6; ++k)
        if (!numbers[k].has_value())
            throw std::invalid_argument("--cutaway expects six numbers, none empty");
    for (std::size_t k = 0; k < 3; ++k) {
        cut.mPoint[k] = *numbers[k];
        cut.mNormal[k] = *numbers[3 + k];
    }
    return cut;
}

RenderOptions cli_render_options(const cli_parsed& rP) {
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
    o.mStreamlines = opt_value(rP, "streamlines");
    if (has_opt(rP, "stream-seeds"))
        o.mStreamSeeds = std::stoi(opt_value(rP, "stream-seeds"));
    if (has_opt(rP, "stream-length"))
        o.mStreamLength = stod_c(opt_value(rP, "stream-length"));
    if (has_opt(rP, "stream-color"))
        o.mStreamColor = cli_parse_rgba(opt_value(rP, "stream-color"), "stream-color");
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

    for (const std::string& text : opt_values(rP, "cutaway"))
        o.mCutaways.push_back(cli_parse_cutaway(text));
    if (has_opt(rP, "cutaway-tint")) {
        if (o.mCutaways.empty())
            throw std::runtime_error("--cutaway-tint requires --cutaway");
        o.mCutawayTint = cli_parse_rgba(opt_value(rP, "cutaway-tint"), "cutaway-tint");
    }

    if (has_opt(rP, "theme")) {
        const std::string theme = opt_value(rP, "theme");
        if (theme == "synthwave")
            o.mTheme = RenderTheme::Synthwave;
        else if (theme != "none")
            throw std::runtime_error("--theme expects synthwave or none, not '" + theme + "'");
    }
    o.mScanlines = has_flag(rP, "scanlines");
    o.mBloom = has_flag(rP, "bloom");
    o.mFringe = has_flag(rP, "fringe");
    o.mGridFloor = has_flag(rP, "grid-floor");
    if (o.mGridFloor && o.mTheme == RenderTheme::None)
        throw std::runtime_error("--grid-floor requires --theme");
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
        // The theme brings its own colormap unless one is named.
        o.mCmap = opt_value(rP, "cmap", o.mTheme == RenderTheme::Synthwave ? "synthwave" : "viridis");
        if (has_opt(rP, "nan-color"))
            o.mNanColor = cli_parse_rgba(opt_value(rP, "nan-color"), "nan-color");
        o.mColorbar = has_flag(rP, "colorbar");
    }
    return o;
}

// ---------------------------------------------------------------------------
// RenderFlags: the flags as a set, and RenderOptions back into flags
// ---------------------------------------------------------------------------

namespace {

// Shortest text that reads back to the same double.
std::string flag_number(double Value) {
    // The shortest of %.15g, %.16g and %.17g that reads back to the same double,
    // so 1.4 stays 1.4 and not 1.3999999999999999.
    char buffer[40];
    for (int digits = 15; digits <= 17; ++digits) {
        detail::snprintf_c(buffer, sizeof(buffer), "%.*g", digits, Value);
        const char* end = nullptr;
        if (detail::parse_double(buffer, end) == Value)
            break;
    }
    return buffer;
}

std::string flag_color(const RenderColor& rColor) {
    static const char digits[] = "0123456789abcdef";
    std::string out = "#";
    for (std::uint8_t c : rColor) {
        out.push_back(digits[c >> 4]);
        out.push_back(digits[c & 15]);
    }
    return out;
}

std::string flag_numbers(const std::vector<double>& rValues) {
    std::string out;
    for (double v : rValues)
        out += (out.empty() ? "" : ",") + flag_number(v);
    return out;
}

const char* edges_text(RenderEdges Edges) {
    return Edges == RenderEdges::All ? "all" : Edges == RenderEdges::Feature ? "feature" : "none";
}

const char* scale_text(RenderScale Scale) {
    return Scale == RenderScale::Log ? "log" : Scale == RenderScale::Symlog ? "symlog" : "linear";
}

const char* diagnostic_text(RenderDiagnostic D) {
    switch (D) {
        case RenderDiagnostic::Quality:
            return "quality";
        case RenderDiagnostic::Inverted:
            return "inverted";
        case RenderDiagnostic::Degenerate:
            return "degenerate";
        case RenderDiagnostic::Orientation:
            return "orientation";
        case RenderDiagnostic::FreeEdges:
            return "free-edges";
        case RenderDiagnostic::EdgeLength:
            return "edge-length";
        default:
            return "none";
    }
}

}  // namespace

RenderFlags RenderFlags::FromOptions(const RenderOptions& rOptions) {
    const RenderOptions d;
    RenderFlags f;
    auto number = [&](const char* pName, double Value, double Default) {
        if (Value != Default)
            f.Set(pName, flag_number(Value));
    };
    auto color = [&](const char* pName, const RenderColor& rValue, const RenderColor& rDefault) {
        if (rValue != rDefault)
            f.Set(pName, flag_color(rValue));
    };
    // Camera. A named view and explicit angles exclude each other; the angles win
    // when both are set, which is what the renderer would draw only if the view
    // were empty, so a view is kept and the angles dropped.
    if (!rOptions.mView.empty()) {
        f.Set("view", rOptions.mView);
    } else {
        number("azimuth", rOptions.mAzimuth, d.mAzimuth);
        number("elevation", rOptions.mElevation, d.mElevation);
    }
    number("roll", rOptions.mRoll, d.mRoll);
    if (rOptions.mProjection == RenderProjection::Perspective) {
        f.SetFlag("perspective", true);
        number("fov", rOptions.mFovDeg, d.mFovDeg);
    }
    number("zoom", rOptions.mZoom, d.mZoom);
    number("pan-x", rOptions.mPanX, d.mPanX);
    number("pan-y", rOptions.mPanY, d.mPanY);
    if (rOptions.mShading != d.mShading)
        f.Set("shading", rOptions.mShading == RenderShading::None ? "none" : "smooth");
    if (!rOptions.mTwoSided)
        f.SetFlag("one-sided", true);
    number("ambient", rOptions.mAmbient, d.mAmbient);
    number("split-angle", rOptions.mSplitAngle, d.mSplitAngle);
    if (rOptions.mEdges != d.mEdges)
        f.Set("edges", edges_text(rOptions.mEdges));
    number("feature-angle", rOptions.mFeatureAngle, d.mFeatureAngle);
    color("edge-color", rOptions.mEdgeColor, d.mEdgeColor);
    color("fill", rOptions.mFillColor, d.mFillColor);
    color("line-color", rOptions.mLineColor, d.mLineColor);
    color("background", rOptions.mBackground, d.mBackground);
    number("point-radius", rOptions.mPointRadius, d.mPointRadius);
    if (rOptions.mSupersample != d.mSupersample)
        f.Set("supersample", std::to_string(rOptions.mSupersample));
    if (rOptions.mAxes)
        f.SetFlag("axes", true);
    if (rOptions.mScaleBar)
        f.SetFlag("scale-bar", true);

    // Field rendering: the range flags only exist beside a mapped field.
    if (!rOptions.mReduce.empty())
        f.Set("reduce", rOptions.mReduce);
    if (!rOptions.mExpr.empty())
        f.Set("expr", rOptions.mExpr);
    if (!rOptions.mColorBy.empty())
        f.Set("color-by", rOptions.mColorBy);
    if (rOptions.mDiagnostic != RenderDiagnostic::None)
        f.Set("diagnostic", diagnostic_text(rOptions.mDiagnostic));
    if (!rOptions.mQualityMetric.empty())
        f.Set("quality-metric", rOptions.mQualityMetric);
    const bool mapped = !rOptions.mColorBy.empty() || !rOptions.mExpr.empty() ||
                        rOptions.mDiagnostic == RenderDiagnostic::Quality ||
                        rOptions.mDiagnostic == RenderDiagnostic::EdgeLength;
    if (mapped) {
        if (rOptions.mComponent.has_value())
            f.Set("component", std::to_string(*rOptions.mComponent));
        if (rOptions.mCmap !=
            (rOptions.mTheme == RenderTheme::Synthwave ? std::string("synthwave") : d.mCmap))
            f.Set("cmap", rOptions.mCmap);
        if (rOptions.mVMin.has_value())
            f.Set("vmin", flag_number(*rOptions.mVMin));
        if (rOptions.mVMax.has_value())
            f.Set("vmax", flag_number(*rOptions.mVMax));
        color("nan-color", rOptions.mNanColor, d.mNanColor);
        if (rOptions.mColorbar)
            f.SetFlag("colorbar", true);
        if (rOptions.mClipLow.has_value() || rOptions.mClipHigh.has_value())
            f.Set("clip", (rOptions.mClipLow ? flag_number(*rOptions.mClipLow) : "") + "," +
                              (rOptions.mClipHigh ? flag_number(*rOptions.mClipHigh) : ""));
        if (rOptions.mSymmetric)
            f.SetFlag("symmetric", true);
        if (rOptions.mScale != d.mScale)
            f.Set("scale", scale_text(rOptions.mScale));
        number("scale-threshold", rOptions.mScaleThreshold, d.mScaleThreshold);
    }
    if (rOptions.mCategorical)
        f.SetFlag("categorical", true);
    if (rOptions.mColorRegions)
        f.SetFlag("color-regions", true);
    if (rOptions.mCategoryEdges)
        f.SetFlag("category-edges", true);
    if (rOptions.mIsolines != d.mIsolines)
        f.Set("isolines", std::to_string(rOptions.mIsolines));
    if (!rOptions.mIsoLevels.empty())
        f.Set("iso-levels", flag_numbers(rOptions.mIsoLevels));
    color("iso-color", rOptions.mIsoColor, d.mIsoColor);
    if (!rOptions.mVectors.empty())
        f.Set("vectors", rOptions.mVectors);
    if (rOptions.mVectorCount != d.mVectorCount)
        f.Set("vector-count", std::to_string(rOptions.mVectorCount));
    number("vector-length", rOptions.mVectorLength, d.mVectorLength);
    color("vector-color", rOptions.mVectorColor, d.mVectorColor);
    if (!rOptions.mStreamlines.empty())
        f.Set("streamlines", rOptions.mStreamlines);
    if (rOptions.mStreamSeeds != d.mStreamSeeds)
        f.Set("stream-seeds", std::to_string(rOptions.mStreamSeeds));
    number("stream-length", rOptions.mStreamLength, d.mStreamLength);
    color("stream-color", rOptions.mStreamColor, d.mStreamColor);
    if (!rOptions.mWarp.empty())
        f.Set("warp", rOptions.mWarp);
    number("warp-scale", rOptions.mWarpScale, d.mWarpScale);
    if (rOptions.mWarpOutline)
        f.SetFlag("warp-outline", true);
    color("outline-color", rOptions.mOutlineColor, d.mOutlineColor);
    for (const RenderCutaway& cut : rOptions.mCutaways)
        f.Add("cutaway", flag_numbers({cut.mPoint[0], cut.mPoint[1], cut.mPoint[2], cut.mNormal[0],
                                       cut.mNormal[1], cut.mNormal[2]}));
    if (!rOptions.mCutaways.empty())
        color("cutaway-tint", rOptions.mCutawayTint, d.mCutawayTint);
    if (rOptions.mTheme == RenderTheme::Synthwave)
        f.Set("theme", "synthwave");
    if (rOptions.mScanlines)
        f.SetFlag("scanlines", true);
    if (rOptions.mBloom)
        f.SetFlag("bloom", true);
    if (rOptions.mFringe)
        f.SetFlag("fringe", true);
    if (rOptions.mGridFloor)
        f.SetFlag("grid-floor", true);
    return f;
}

RenderFlags RenderFlags::FromTokens(const std::vector<std::string>& rTokens) {
    const cli_parsed parsed = cli_parse(rTokens, render_flag_specs());
    RenderFlags f;
    for (const auto& [name, value] : parsed.ordered)
        f.Add(name, value);
    for (const std::string& name : parsed.flags)
        f.mFlags.insert(name);
    // A repeated value flag other than `cutaway` keeps its last value, as the
    // parser does.
    for (auto& [name, values] : f.mValues)
        if (name != "cutaway" && values.size() > 1)
            values = {values.back()};
    return f;
}

bool RenderFlags::Known(const std::string& rName) {
    for (const cli_opt_spec& spec : render_flag_specs())
        if (spec.canonical == rName)
            return true;
    return false;
}

bool RenderFlags::TakesValue(const std::string& rName) {
    for (const cli_opt_spec& spec : render_flag_specs())
        if (spec.canonical == rName)
            return spec.takes_value;
    return false;
}

bool RenderFlags::Has(const std::string& rName) const {
    return mValues.count(rName) != 0 || mFlags.count(rName) != 0;
}

std::string RenderFlags::Value(const std::string& rName, const std::string& rDefault) const {
    const auto it = mValues.find(rName);
    return it == mValues.end() || it->second.empty() ? rDefault : it->second.back();
}

const std::vector<std::string>& RenderFlags::Values(const std::string& rName) const {
    static const std::vector<std::string> empty;
    const auto it = mValues.find(rName);
    return it == mValues.end() ? empty : it->second;
}

void RenderFlags::Set(const std::string& rName, const std::string& rValue) {
    mValues[rName] = {rValue};
}

void RenderFlags::Add(const std::string& rName, const std::string& rValue) {
    mValues[rName].push_back(rValue);
}

void RenderFlags::SetFlag(const std::string& rName, bool On) {
    if (On)
        mFlags.insert(rName);
    else
        mFlags.erase(rName);
}

void RenderFlags::Unset(const std::string& rName) {
    mValues.erase(rName);
    mFlags.erase(rName);
}

RenderOptions RenderFlags::ToOptions() const {
    cli_parsed parsed;
    for (const auto& [name, values] : mValues)
        for (const std::string& value : values) {
            parsed.values[name] = value;
            parsed.multi[name].push_back(value);
            parsed.ordered.emplace_back(name, value);
        }
    parsed.flags.insert(mFlags.begin(), mFlags.end());
    return cli_render_options(parsed);
}

std::vector<std::string> RenderFlags::Tokens() const {
    std::vector<std::string> names;
    for (const auto& entry : mValues)
        names.push_back(entry.first);
    for (const std::string& name : mFlags)
        names.push_back(name);
    std::sort(names.begin(), names.end());
    std::vector<std::string> out;
    for (const std::string& name : names) {
        const auto it = mValues.find(name);
        if (it != mValues.end()) {
            for (const std::string& value : it->second) {
                // `--name=value`: a value that starts with a dash survives parsing.
                out.push_back("--" + name + "=" + value);
            }
        } else {
            out.push_back("--" + name);
        }
    }
    return out;
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
