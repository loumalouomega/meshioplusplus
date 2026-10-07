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
 * @file operations/render.hpp
 * @brief A deterministic software rasterizer, and the terminal, HTML and PNG
 * encodings of its frames.
 *
 * `render` turns a mesh into an RGBA frame with a depth buffer, so faces that
 * intersect or overlap are resolved per pixel -- the step the painter's
 * algorithm of the SVG and TikZ writers cannot take. The frame then becomes
 * text (`encode_text`: Unicode block, sextant, Braille or ASCII cells with
 * 24-bit, 256-, 16-colour or no SGR colour, as plain text, ANSI or an HTML
 * `<pre>`), a graphics-protocol escape sequence (Kitty, iTerm2, Sixel), or a
 * PNG (`encode_png`, dependency-free). `write_snapshot` does all of it from a
 * path's extension. See `doc/tui.md`.
 *
 * ### What is drawn
 *
 * The faces the SVG writer draws: 2-D cells as they are, volume cells (and
 * polyhedra) through their boundary skin, `line` cells as lines and `vertex`
 * cells -- or the points of a mesh with no cells -- as discs. Higher-order
 * cells are drawn through their corners. A volume mesh's skin and its own 2-D
 * cells are both drawn; where they coincide the earlier primitive (the skin)
 * wins the depth tie.
 *
 * ### Determinism
 *
 * A frame is a pure function of the mesh and the options: the same bytes on
 * every parallel backend and thread count. Screen positions are fixed point
 * with 8 subpixel bits, coverage is decided by integer edge functions with a
 * top-left fill rule, the depth test is strict so the earlier primitive (in
 * block-major order) wins a tie, each 32x32 tile is owned by one worker, and
 * supersampling resolves through a fixed-order integer box filter.
 *
 * Colouring follows `detail/face_color.hpp`, the resolver the SVG and TikZ
 * writers share: point data colours a face by the mean of its corners, cell
 * data by the owning cell, and the automatic range spans the drawn values.
 */

// System includes
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The camera model.
enum class RenderProjection : std::uint8_t {
    Orthographic = 0,  ///< parallel rays; the SVG writer's camera
    Perspective = 1,   ///< a pinhole with `mFovDeg` vertical field of view
};

/// How faces are lit. Shading multiplies the field (or fill) colour; it never
/// replaces it.
enum class RenderShading : std::uint8_t {
    None = 0,    ///< unlit: every face shows its colour as is
    Flat = 1,    ///< one Lambert term per face from its Newell normal (default)
    Smooth = 2,  ///< per-pixel normals interpolated from `compute_normals`, split at `mSplitAngle`
};

/// Which edges are drawn over the faces.
enum class RenderEdges : std::uint8_t {
    None = 0,     ///< no overlay (default)
    All = 1,      ///< every element edge (polygon sides, not triangulation diagonals)
    Feature = 2,  ///< the sharp, open and non-manifold edges of `feature_edges`
};

/// An RGBA colour, 8 bits per channel, alpha 255 opaque.
using RenderColor = std::array<std::uint8_t, 4>;

/** @brief Everything `render` takes besides the mesh. */
struct RenderOptions {
    /// Frame size in pixels; both must be positive.
    int mWidth = 320;
    int mHeight = 240;
    /// Height over width of one pixel: 1 for an image, other values for the
    /// non-square subpixels of a terminal cell (`text_frame_size` sets it).
    double mPixelAspect = 1.0;
    /// Supersampling factor per axis: 1, 2 or 4.
    int mSupersample = 1;

    /// Camera direction in degrees, with `camera_basis`'s conventions (the SVG
    /// writer's): the defaults are the isometric view.
    double mAzimuth = 45.0;
    double mElevation = 35.264389682754654;
    double mRoll = 0.0;
    /// A named view that overrides azimuth and elevation when non-empty:
    /// `iso`, `+x`, `-x`, `+y`, `-y`, `+z`, `-z` (the camera sits on that side).
    std::string mView;
    RenderProjection mProjection = RenderProjection::Orthographic;
    /// Vertical field of view of the perspective camera, in (0, 180).
    double mFovDeg = 30.0;
    /// Magnification over the fit-to-frame view; must be positive.
    double mZoom = 1.0;
    /// Shift of the view, in fractions of the frame width and height.
    double mPanX = 0.0;
    double mPanY = 0.0;

    RenderShading mShading = RenderShading::Flat;
    /// Light both sides of a face, so an open shell's inside and an inverted
    /// element are visible; one-sided lighting leaves back faces at `mAmbient`.
    bool mTwoSided = true;
    /// The light floor in [0, 1]: a face turned away from the light keeps this
    /// fraction of its colour.
    double mAmbient = 0.25;
    /// Direction toward the light in camera space (x right, y up, z toward the
    /// viewer); the default is a headlight.
    std::array<double, 3> mLightDir = {0.0, 0.0, 1.0};
    /// The crease angle for `Smooth` shading, in [0, 180] degrees.
    double mSplitAngle = 30.0;

    RenderEdges mEdges = RenderEdges::None;
    /// The feature angle for `RenderEdges::Feature`, in [0, 180] degrees.
    double mFeatureAngle = 30.0;
    RenderColor mEdgeColor = {0, 0, 0, 255};
    /// Colour of faces when no field is mapped (the SVG writer's `#c8c5bd`).
    RenderColor mFillColor = {200, 197, 189, 255};
    /// Colour of `line` cells when no field is mapped (the SVG stroke, `#000080`).
    RenderColor mLineColor = {0, 0, 128, 255};
    /// The background; alpha 0 (the default) is transparent in a PNG and the
    /// terminal's own background in text.
    RenderColor mBackground = {0, 0, 0, 0};
    /// Disc radius of drawn points, in output pixels.
    double mPointRadius = 1.5;

    /// Data array to colour by; empty for the fill colour.
    std::string mColorBy;
    /// Component of a multi-component array; its magnitude when unset.
    std::optional<int> mComponent;
    /// Built-in colormap name (`colormap_names()`).
    std::string mCmap = "viridis";
    std::optional<double> mVMin;
    std::optional<double> mVMax;
    /// Colour of faces whose value is NaN or infinite.
    RenderColor mNanColor = {128, 128, 128, 255};
    /// Draw a colour bar along the right edge and name the range in `Frame::mNotes`.
    bool mColorbar = false;

    /// Draw the world axes (x red, y green, z blue) in the lower-left corner.
    bool mAxes = false;
    /// Draw a scale bar in the lower-right corner (orthographic only) and give
    /// its length in `Frame::mNotes`.
    bool mScaleBar = false;
};

/** @brief A rendered image. */
struct Frame {
    int mWidth = 0;
    int mHeight = 0;
    /// `mWidth * mHeight * 4` bytes, rows top to bottom, straight (not
    /// premultiplied) alpha.
    std::vector<std::uint8_t> mRgba;
    /// Per pixel, the input cell drawn there (global block-major index, the
    /// convention of regions and `"surface:parent_cell"`); -1 for the
    /// background, overlays and the points of a cell-less cloud.
    std::vector<std::int64_t> mCellIds;
    /// Whether a field was mapped, and the range it was mapped over.
    bool mColored = false;
    double mVMin = 0.0;
    double mVMax = 0.0;
    /// Lines of text describing the frame (the colour range, the scale bar),
    /// printed under a text rendering.
    std::vector<std::string> mNotes;
};

/// How a frame's pixels become terminal cells, or a graphics protocol.
enum class TextEncoding : std::uint8_t {
    HalfBlock = 0,  ///< 1x2 pixels per cell, `▀`/`▄` (default)
    Quadrant = 1,   ///< 2x2, the quadrant block elements
    Sextant = 2,    ///< 2x3, Unicode 13 block sextants
    Braille = 3,    ///< 2x4 Braille dots
    Ascii = 4,      ///< 1x1, a luminance ramp of printable ASCII
    Kitty = 5,      ///< the Kitty graphics protocol (RGBA, chunked base64)
    ITerm2 = 6,     ///< iTerm2 inline images (OSC 1337, a PNG)
    Sixel = 7,      ///< DEC Sixel, a 256-colour median-cut palette
};

/// How many colours a cell encoding may use.
enum class ColorDepth : std::uint8_t {
    TrueColor = 0,   ///< 24-bit SGR `38;2;r;g;b`
    Palette256 = 1,  ///< the xterm 256-colour palette
    Ansi16 = 2,      ///< the 16 ANSI colours
    Mono = 3,        ///< no colour: glyphs only
};

/// The container of a cell encoding.
enum class TextFormat : std::uint8_t {
    Ansi = 0,   ///< glyphs with SGR colour escapes, for a terminal or `cat`
    Plain = 1,  ///< glyphs only, no escape sequence at all
    Html = 2,   ///< a self-contained HTML page holding a `<pre>` of coloured spans
};

/** @brief Options of `encode_text` and `render_text`. */
struct TextOptions {
    TextEncoding mEncoding = TextEncoding::HalfBlock;
    ColorDepth mDepth = ColorDepth::TrueColor;
    TextFormat mFormat = TextFormat::Ansi;
    /// Terminal size `render_text` fits the frame to, in cells.
    int mCols = 80;
    int mRows = 24;
    /// Height over width of one terminal cell.
    double mCellAspect = 2.0;
    /// Pixel size of one cell, used to size a graphics-protocol image.
    int mCellPixelWidth = 8;
    int mCellPixelHeight = 16;
    /// Wrap graphics-protocol sequences for tmux's passthrough
    /// (`allow-passthrough on`, tmux 3.3 and later).
    bool mTmuxPassthrough = false;
    /// Append `Frame::mNotes` under the picture.
    bool mNotes = true;
};

/** @brief Render a mesh into an RGBA frame. */
MESHIOPLUSPLUS_API Frame render(const Mesh& rMesh, const RenderOptions& rOptions = {});

/**
 * @brief The frame size, in pixels, and the pixel aspect that make an encoding
 * fill `Cols` by `Rows` cells.
 * @return `{width, height}`; @p rPixelAspect receives the matching
 *         `RenderOptions::mPixelAspect`
 */
MESHIOPLUSPLUS_API std::array<int, 2> text_frame_size(const TextOptions& rOptions,
                                                      double& rPixelAspect);

/**
 * @brief Encode a frame as text.
 * A cell encoding needs a frame whose size is a multiple of its cell grid
 * (`text_frame_size` gives one); each cell takes the two-colour split of its
 * pixels with the least squared error, ties broken by a fixed order. A
 * graphics protocol takes any frame size and needs `TextFormat::Ansi`.
 * @throws std::invalid_argument on a frame the encoding cannot take
 */
MESHIOPLUSPLUS_API std::string encode_text(const Frame& rFrame, const TextOptions& rOptions = {});

/** @brief `render` sized by `text_frame_size`, then `encode_text`. */
MESHIOPLUSPLUS_API std::string render_text(const Mesh& rMesh, const RenderOptions& rRender = {},
                                           const TextOptions& rText = {});

/**
 * @brief Encode a frame as an RGBA PNG.
 * @param CompressLevel 0 (the default) writes stored deflate blocks with no
 *        dependency, the same bytes on every platform; 1-9 compresses through
 *        zlib, whose output can differ between zlib versions.
 * @throws std::invalid_argument for a level outside [0, 9], or 1-9 in a build
 *         without zlib
 */
MESHIOPLUSPLUS_API std::string encode_png(const Frame& rFrame, int CompressLevel = 0);

/**
 * @brief The colour depth a terminal advertises through its environment:
 * `NO_COLOR` set (to anything) gives `Mono`; `COLORTERM` of `truecolor` or
 * `24bit` gives `TrueColor`; a `TERM` containing `256color` gives
 * `Palette256`; `TERM` of `dumb` gives `Mono`; anything else `Ansi16`.
 * Null pointers mean unset variables.
 */
MESHIOPLUSPLUS_API ColorDepth detect_color_depth(const char* pNoColor, const char* pColorTerm,
                                                 const char* pTerm);

/** @brief Options of `write_snapshot` beyond the render and text ones. */
struct SnapshotOptions {
    /// PNG compression level (see `encode_png`).
    int mPngCompress = 0;
    /// An asciicast orbit: frames, rate, and the azimuth swept.
    int mCastFrames = 36;
    double mCastFps = 12.0;
    double mCastDegrees = 360.0;
};

/**
 * @brief The asciicast v2 recording of an orbit around the mesh: frame `i` is
 * rendered at azimuth `mAzimuth + i * mCastDegrees / mCastFrames` and stamped
 * `i / mCastFps` seconds, never the wall clock, so the file is reproducible.
 */
MESHIOPLUSPLUS_API std::string encode_cast(const Mesh& rMesh, const RenderOptions& rRender,
                                           const TextOptions& rText,
                                           const SnapshotOptions& rSnapshot = {});

/**
 * @brief Render a mesh to a file chosen by extension: `.png` (`encode_png` of
 * `render`), `.txt` (plain cells), `.ansi` (ANSI cells), `.html` (an HTML
 * page) or `.cast` (`encode_cast`). For the text forms the frame is sized by
 * `rText`'s columns and rows; `rText.mFormat` is overridden by the extension.
 * @throws std::invalid_argument for any other extension
 */
MESHIOPLUSPLUS_API void write_snapshot(const std::string& rPath, const Mesh& rMesh,
                                       const RenderOptions& rRender = {},
                                       const TextOptions& rText = {},
                                       const SnapshotOptions& rSnapshot = {});

}  // namespace meshioplusplus
