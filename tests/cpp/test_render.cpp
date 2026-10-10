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

// The software rasterizer and its encodings (operations/render.hpp): the
// probes of roadmap 7.2.7 and 7.3.6, plus the PNG and asciicast containers.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/operations/render.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/png_write.hpp"
#include "../../src/cpp/src/detail/raster.hpp"
#include "../../src/cpp/src/detail/render_field.hpp"
#include "../../src/cpp/src/detail/streamlines.hpp"

using namespace meshioplusplus;

namespace {

using mt::make_mesh;

Mesh cube_quads() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}});
}

std::set<std::int64_t> ids_of(const Frame& rFrame) {
    return std::set<std::int64_t>(rFrame.mCellIds.begin(), rFrame.mCellIds.end());
}

Frame solid_frame(int W, int H, std::uint8_t R, std::uint8_t G, std::uint8_t B,
                  std::uint8_t A = 255) {
    Frame f;
    f.mWidth = W;
    f.mHeight = H;
    for (int i = 0; i < W * H; ++i) {
        f.mRgba.push_back(R);
        f.mRgba.push_back(G);
        f.mRgba.push_back(B);
        f.mRgba.push_back(A);
    }
    f.mCellIds.assign(static_cast<std::size_t>(W * H), -1);
    return f;
}

void set_px(Frame& rFrame, int X, int Y, std::uint8_t R, std::uint8_t G, std::uint8_t B,
            std::uint8_t A = 255) {
    std::uint8_t* p = rFrame.mRgba.data() + (static_cast<std::size_t>(Y) * rFrame.mWidth + X) * 4;
    p[0] = R;
    p[1] = G;
    p[2] = B;
    p[3] = A;
}

TextOptions text(TextEncoding Encoding, ColorDepth Depth, TextFormat Format = TextFormat::Ansi) {
    TextOptions t;
    t.mEncoding = Encoding;
    t.mDepth = Depth;
    t.mFormat = Format;
    return t;
}

}  // namespace

// ---------------------------------------------------------------------------
// The raster core (7.2.2, 7.2.7)
// ---------------------------------------------------------------------------

// Corners on pixel centres: the top edge (y = 0.5, running right) and the left
// edge (x = 0.5, running up) own the centres on them; the hypotenuse
// x + y = 4 does not. So exactly the six pixels with x + y <= 2 are covered,
// whichever way round the corners are given.
TEST(Render, TopLeftRuleCoversTheHandDerivedPixels) {
    for (bool flip : {false, true}) {
        detail::RasterScene scene;
        scene.mVertices = {{0.5, 0.5, 0.0}, {3.5, 0.5, 0.0}, {0.5, 3.5, 0.0}};
        detail::RasterTri tri;
        tri.mV[0] = 0;
        tri.mV[1] = flip ? 2 : 1;
        tri.mV[2] = flip ? 1 : 2;
        tri.mId = 7;
        scene.mTris.push_back(tri);
        const detail::RasterTarget t = detail::rasterize(scene, 4, 4, {0, 0, 0, 0});
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
                EXPECT_EQ(t.mIds[static_cast<std::size_t>(y * 4 + x)], x + y <= 2 ? 7 : -1)
                    << "pixel " << x << "," << y << " flip " << flip;
    }
}

// Two triangles sharing a diagonal through pixel centres cover every pixel
// exactly once: no gap and no pixel claimed twice (the second would lose the
// strict depth test, so a double claim shows as the first id everywhere).
TEST(Render, SharedEdgeIsOwnedByExactlyOneTriangle) {
    detail::RasterScene scene;
    scene.mVertices = {{0.0, 0.0, 0.0}, {8.0, 0.0, 0.0}, {8.0, 8.0, 0.0}, {0.0, 8.0, 0.0}};
    detail::RasterTri a;
    a.mV[0] = 0;
    a.mV[1] = 1;
    a.mV[2] = 2;
    a.mId = 1;
    detail::RasterTri b;
    b.mV[0] = 0;
    b.mV[1] = 2;
    b.mV[2] = 3;
    b.mId = 2;
    scene.mTris = {a, b};
    const detail::RasterTarget t = detail::rasterize(scene, 8, 8, {0, 0, 0, 0});
    int count[3] = {0, 0, 0};
    for (std::int64_t id : t.mIds) {
        ASSERT_GE(id, 1);
        ++count[id];
    }
    // The diagonal's 8 centres go to one side: 36 + 28 either way round.
    EXPECT_EQ(count[1] + count[2], 64);
    EXPECT_TRUE((count[1] == 36 && count[2] == 28) || (count[1] == 28 && count[2] == 36));
}

// Two planes crossing along x = 0: seen from +z, z = x is in front for x > 0
// and z = -x for x < 0. The painter's algorithm draws one of them whole.
TEST(Render, InterpenetratingFacesShowTheirIntersection) {
    Mesh m = make_mesh({{-1, -1, -1},
                        {1, -1, 1},
                        {1, 1, 1},
                        {-1, 1, -1},
                        {-1, -1, 1},
                        {1, -1, -1},
                        {1, 1, -1},
                        {-1, 1, 1}},
                       "quad", {{0, 1, 2, 3}, {4, 5, 6, 7}});
    RenderOptions o;
    o.mWidth = 40;
    o.mHeight = 40;
    o.mView = "+z";
    const Frame f = render(m, o);
    int left_a = 0, left_b = 0, right_a = 0, right_b = 0;
    for (int y = 0; y < 40; ++y)
        for (int x = 0; x < 40; ++x) {
            const std::int64_t id = f.mCellIds[static_cast<std::size_t>(y * 40 + x)];
            if (x < 18) {
                left_a += id == 0;
                left_b += id == 1;
            } else if (x >= 22) {
                right_a += id == 0;
                right_b += id == 1;
            }
        }
    EXPECT_GT(right_a, 0);
    EXPECT_EQ(right_b, 0);
    EXPECT_GT(left_b, 0);
    EXPECT_EQ(left_a, 0);
}

// From the isometric camera (on the +x +y +z side) only the faces x = 1,
// y = 1 and z = 1 can own a pixel.
TEST(Render, HiddenFacesOfACubeNeverOwnAPixel) {
    RenderOptions o;
    o.mWidth = 64;
    o.mHeight = 64;
    for (int ss : {1, 2, 4}) {
        o.mSupersample = ss;
        const std::set<std::int64_t> ids = ids_of(render(cube_quads(), o));
        EXPECT_EQ(ids, (std::set<std::int64_t>{-1, 1, 3, 5})) << "supersample " << ss;
    }
}

TEST(Render, VolumeCellsDrawTheirSkinAndReportTheParentCell) {
    RenderOptions o;
    o.mWidth = 32;
    o.mHeight = 32;
    const std::set<std::int64_t> ids = ids_of(render(mt::hex_mesh(), o));
    EXPECT_EQ(ids, (std::set<std::int64_t>{-1, 0}));
}

TEST(Render, DegenerateInputsKeepTheFrameSize) {
    RenderOptions o;
    o.mWidth = 13;
    o.mHeight = 7;
    Mesh empty;
    for (const Mesh* p : {&empty}) {
        const Frame f = render(*p, o);
        EXPECT_EQ(f.mWidth, 13);
        EXPECT_EQ(f.mHeight, 7);
        EXPECT_EQ(f.mRgba.size(), 13u * 7u * 4u);
        EXPECT_EQ(ids_of(f), (std::set<std::int64_t>{-1}));
    }
    Mesh point;
    point.AssignPoints(mt::points_from({{1.0, 2.0, 3.0}}));
    const Frame fp = render(point, o);
    EXPECT_EQ(fp.mRgba.size(), 13u * 7u * 4u);
    bool drawn = false;
    for (std::size_t i = 3; i < fp.mRgba.size(); i += 4)
        drawn = drawn || fp.mRgba[i] != 0;
    EXPECT_TRUE(drawn);
    const Mesh one = make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}});
    EXPECT_EQ(render(one, o).mRgba.size(), 13u * 7u * 4u);
}

TEST(Render, IsAPureFunctionOfItsInput) {
    RenderOptions o;
    o.mWidth = 50;
    o.mHeight = 30;
    o.mSupersample = 2;
    o.mEdges = RenderEdges::All;
    o.mShading = RenderShading::Smooth;
    const Frame a = render(mt::hex_mesh(), o);
    const Frame b = render(mt::hex_mesh(), o);
    EXPECT_EQ(a.mRgba, b.mRgba);
    EXPECT_EQ(a.mCellIds, b.mCellIds);
}

TEST(Render, ColoursByCellDataAndReportsTheRange) {
    Mesh m = cube_quads();
    m.AddCellData("v", {mt::points_from({{0}, {1}, {2}, {3}, {4}, {5}})});
    RenderOptions o;
    o.mWidth = 32;
    o.mHeight = 32;
    o.mColorBy = "v";
    o.mShading = RenderShading::None;
    const Frame f = render(m, o);
    EXPECT_TRUE(f.mColored);
    // The auto range spans every face handed to the rasterizer, hidden or not
    // (the SVG writer's "drawn faces"), so an orbit never changes the colours.
    EXPECT_DOUBLE_EQ(f.mVMin, 0.0);
    EXPECT_DOUBLE_EQ(f.mVMax, 5.0);
    ASSERT_EQ(f.mNotes.size(), 1u);
    EXPECT_EQ(f.mNotes[0], "v: 0 .. 5 (viridis)");
    o.mColorBy = "nope";
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(Render, ConstantFieldGivesOneColour) {
    Mesh m = cube_quads();
    m.AddPointData("c", mt::points_from({{2}, {2}, {2}, {2}, {2}, {2}, {2}, {2}}));
    RenderOptions o;
    o.mWidth = 24;
    o.mHeight = 24;
    o.mColorBy = "c";
    o.mShading = RenderShading::None;
    const Frame f = render(m, o);
    std::set<std::uint32_t> colors;
    for (std::size_t i = 0; i < f.mCellIds.size(); ++i)
        if (f.mCellIds[i] >= 0)
            colors.insert((std::uint32_t(f.mRgba[i * 4]) << 16) |
                          (std::uint32_t(f.mRgba[i * 4 + 1]) << 8) | f.mRgba[i * 4 + 2]);
    EXPECT_EQ(colors.size(), 1u);
}

TEST(Render, RejectsBadOptions) {
    RenderOptions o;
    o.mSupersample = 3;
    EXPECT_THROW(render(cube_quads(), o), std::invalid_argument);
    o = RenderOptions();
    o.mView = "sideways";
    EXPECT_THROW(render(cube_quads(), o), std::invalid_argument);
    o = RenderOptions();
    o.mWidth = 0;
    EXPECT_THROW(render(cube_quads(), o), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Cell encodings (7.3.1, 7.3.2, 7.3.6)
// ---------------------------------------------------------------------------

TEST(RenderText, SolidFramesEncodeToTheirColourAtEveryDepth) {
    const Frame black = solid_frame(1, 2, 0, 0, 0);
    const Frame white = solid_frame(1, 2, 255, 255, 255);
    const TextEncoding hb = TextEncoding::HalfBlock;
    EXPECT_EQ(encode_text(black, text(hb, ColorDepth::TrueColor)), "\x1b[48;2;0;0;0m \x1b[0m\n");
    EXPECT_EQ(encode_text(white, text(hb, ColorDepth::TrueColor)),
              "\x1b[48;2;255;255;255m \x1b[0m\n");
    EXPECT_EQ(encode_text(black, text(hb, ColorDepth::Palette256)), "\x1b[48;5;16m \x1b[0m\n");
    EXPECT_EQ(encode_text(white, text(hb, ColorDepth::Palette256)), "\x1b[48;5;231m \x1b[0m\n");
    EXPECT_EQ(encode_text(black, text(hb, ColorDepth::Ansi16)), "\x1b[40m \x1b[0m\n");
    EXPECT_EQ(encode_text(white, text(hb, ColorDepth::Ansi16)), "\x1b[107m \x1b[0m\n");
    EXPECT_EQ(encode_text(black, text(hb, ColorDepth::Mono)), " \n");
    EXPECT_EQ(encode_text(white, text(hb, ColorDepth::Mono)), "\xe2\x96\x88\n");
}

TEST(RenderText, TwoColourCheckerRoundTripsThroughHalfBlocks) {
    Frame f = solid_frame(2, 2, 0, 0, 0);
    set_px(f, 0, 0, 255, 0, 0);
    set_px(f, 1, 0, 0, 0, 255);
    set_px(f, 0, 1, 0, 0, 255);
    set_px(f, 1, 1, 255, 0, 0);
    EXPECT_EQ(encode_text(f, text(TextEncoding::HalfBlock, ColorDepth::TrueColor)),
              "\x1b[38;2;255;0;0;48;2;0;0;255m\xe2\x96\x80"
              "\x1b[38;2;0;0;255;48;2;255;0;0m\xe2\x96\x80\x1b[0m\n");
}

TEST(RenderText, TransparentPixelsUseTheTerminalBackground) {
    Frame f = solid_frame(2, 2, 0, 0, 0, 0);
    set_px(f, 1, 1, 255, 0, 0);
    // Only the lower-right quadrant is drawn: the glyph is drawn in red on the
    // default background, never the transparent group as the foreground.
    EXPECT_EQ(encode_text(f, text(TextEncoding::Quadrant, ColorDepth::TrueColor)),
              "\x1b[38;2;255;0;0m\xe2\x96\x97\x1b[0m\n");
    EXPECT_EQ(encode_text(solid_frame(1, 2, 0, 0, 0, 0),
                          text(TextEncoding::HalfBlock, ColorDepth::TrueColor)),
              " \n");
}

TEST(RenderText, SextantAndBrailleCodePoints) {
    // Only subpixels 1 and 3 (top-left and middle-left... index 0 and 2) lit:
    // mask 0b000101 = 5 -> U+1FB04 (BLOCK SEXTANT-13).
    Frame s = solid_frame(2, 3, 0, 0, 0, 0);
    set_px(s, 0, 0, 255, 255, 255);
    set_px(s, 0, 1, 255, 255, 255);
    EXPECT_EQ(encode_text(s, text(TextEncoding::Sextant, ColorDepth::Mono)), "\xf0\x9f\xac\x84\n");
    // The left column of a Braille cell is dots 1, 2, 3, 7: U+2847.
    Frame b = solid_frame(2, 4, 0, 0, 0, 0);
    for (int y = 0; y < 4; ++y)
        set_px(b, 0, y, 255, 255, 255);
    EXPECT_EQ(encode_text(b, text(TextEncoding::Braille, ColorDepth::Mono)), "\xe2\xa1\x87\n");
}

TEST(RenderText, AsciiRampIsMonotoneInLuminance) {
    const std::string ramp = " .:-=+*#%@";
    std::size_t last = 0;
    for (int v = 0; v < 256; ++v) {
        const std::string out =
            encode_text(solid_frame(1, 1, v, v, v), text(TextEncoding::Ascii, ColorDepth::Mono));
        ASSERT_EQ(out.size(), 2u);
        const std::size_t pos = ramp.find(out[0]);
        ASSERT_NE(pos, std::string::npos);
        EXPECT_GE(pos, last);
        last = pos;
    }
    EXPECT_EQ(last, ramp.size() - 1);
}

TEST(RenderText, NoColorEmitsNoEscapeSequence) {
    EXPECT_EQ(detect_color_depth("1", "truecolor", "xterm-256color"), ColorDepth::Mono);
    EXPECT_EQ(detect_color_depth("", "truecolor", nullptr), ColorDepth::TrueColor);
    EXPECT_EQ(detect_color_depth(nullptr, "24bit", nullptr), ColorDepth::TrueColor);
    EXPECT_EQ(detect_color_depth(nullptr, nullptr, "xterm-256color"), ColorDepth::Palette256);
    EXPECT_EQ(detect_color_depth(nullptr, nullptr, "dumb"), ColorDepth::Mono);
    EXPECT_EQ(detect_color_depth(nullptr, nullptr, "xterm"), ColorDepth::Ansi16);
    RenderOptions o;
    TextOptions t;
    t.mCols = 20;
    t.mRows = 10;
    t.mDepth = detect_color_depth("1", nullptr, nullptr);
    o.mColorBy = "";
    const std::string out = render_text(cube_quads(), o, t);
    EXPECT_EQ(out.find('\x1b'), std::string::npos);
}

TEST(RenderText, FrameSizeFollowsTheEncoding) {
    TextOptions t;
    t.mCols = 10;
    t.mRows = 4;
    double aspect = 0.0;
    t.mEncoding = TextEncoding::HalfBlock;
    EXPECT_EQ(text_frame_size(t, aspect), (std::array<int, 2>{10, 8}));
    EXPECT_DOUBLE_EQ(aspect, 1.0);
    t.mEncoding = TextEncoding::Braille;
    EXPECT_EQ(text_frame_size(t, aspect), (std::array<int, 2>{20, 16}));
    EXPECT_DOUBLE_EQ(aspect, 1.0);
    t.mEncoding = TextEncoding::Sextant;
    EXPECT_EQ(text_frame_size(t, aspect), (std::array<int, 2>{20, 12}));
    EXPECT_NEAR(aspect, 4.0 / 3.0, 1e-15);
    EXPECT_THROW(
        encode_text(solid_frame(3, 3, 0, 0, 0), text(TextEncoding::Quadrant, ColorDepth::Mono)),
        std::invalid_argument);
}

TEST(RenderText, HtmlIsASelfContainedPage) {
    Frame f = solid_frame(1, 2, 255, 0, 0);
    const std::string html =
        encode_text(f, text(TextEncoding::HalfBlock, ColorDepth::TrueColor, TextFormat::Html));
    EXPECT_EQ(html.rfind("<!DOCTYPE html>", 0), 0u);
    EXPECT_NE(html.find("<span style=\"background:#ff0000\"> </span>"), std::string::npos);
    EXPECT_NE(html.find("</html>"), std::string::npos);
}

TEST(RenderText, GraphicsProtocols) {
    Frame f = solid_frame(3, 2, 10, 20, 30);
    const std::string kitty = encode_text(f, text(TextEncoding::Kitty, ColorDepth::TrueColor));
    EXPECT_EQ(kitty.rfind("\x1b_Ga=T,f=32,s=3,v=2,q=2,m=0;", 0), 0u);
    TextOptions tmux = text(TextEncoding::Kitty, ColorDepth::TrueColor);
    tmux.mTmuxPassthrough = true;
    EXPECT_EQ(encode_text(f, tmux).rfind("\x1bPtmux;\x1b\x1b_G", 0), 0u);
    const std::string iterm = encode_text(f, text(TextEncoding::ITerm2, ColorDepth::TrueColor));
    EXPECT_EQ(iterm.rfind("\x1b]1337;File=inline=1;", 0), 0u);
    const std::string sixel = encode_text(f, text(TextEncoding::Sixel, ColorDepth::TrueColor));
    // One colour, 10/20/30 out of 255 as percentages; a 3-wide band of rows 0-1.
    EXPECT_EQ(sixel, "\x1bP0;1;0q\"1;1;3;2#0;2;4;8;12#0BBB-\x1b\\\n");
    EXPECT_THROW(encode_text(f, text(TextEncoding::Sixel, ColorDepth::TrueColor, TextFormat::Html)),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// PNG and the other containers (7.4)
// ---------------------------------------------------------------------------

TEST(RenderPng, ChecksumsMatchTheirReferenceValues) {
    const std::string check = "123456789";
    EXPECT_EQ(
        detail::png_crc32(0, reinterpret_cast<const unsigned char*>(check.data()), check.size()),
        0xCBF43926U);
    const std::string wiki = "Wikipedia";
    EXPECT_EQ(
        detail::png_adler32(1, reinterpret_cast<const unsigned char*>(wiki.data()), wiki.size()),
        0x11E60398U);
}

TEST(RenderPng, StoredBlocksInflateBackToTheRows) {
    Frame f = solid_frame(300, 120, 1, 2, 3, 4);  // > 65535 bytes: several stored blocks
    set_px(f, 299, 119, 9, 8, 7, 6);
    const std::string png = encode_png(f);
    ASSERT_EQ(png.substr(0, 8), std::string("\x89PNG\r\n\x1a\n", 8));
    EXPECT_EQ(png.substr(png.size() - 12), std::string("\0\0\0\0IEND\xae\x42\x60\x82", 12));
    // IHDR: 300 x 120, 8-bit RGBA.
    EXPECT_EQ(png.substr(12, 4), "IHDR");
    EXPECT_EQ(static_cast<unsigned char>(png[19]), 44);  // 300 = 0x012C
    EXPECT_EQ(static_cast<unsigned char>(png[25]), 6);
    if (!detail::zlib_available())
        GTEST_SKIP() << "no zlib to inflate with";
    const std::size_t idat_len =
        (static_cast<unsigned char>(png[33]) << 24) | (static_cast<unsigned char>(png[34]) << 16) |
        (static_cast<unsigned char>(png[35]) << 8) | static_cast<unsigned char>(png[36]);
    ASSERT_EQ(png.substr(37, 4), "IDAT");
    const std::string raw =
        detail::zlib_inflate(std::string_view(png).substr(41, idat_len), 15, nullptr, "test");
    ASSERT_EQ(raw.size(), (300u * 4u + 1u) * 120u);
    EXPECT_EQ(raw[0], 0);
    EXPECT_EQ(static_cast<unsigned char>(raw[raw.size() - 4]), 9);
    EXPECT_EQ(static_cast<unsigned char>(raw[raw.size() - 1]), 6);
}

TEST(RenderPng, IsByteIdenticalAndRejectsBadLevels) {
    const Frame f = render(cube_quads(), RenderOptions());
    EXPECT_EQ(encode_png(f), encode_png(f));
    EXPECT_THROW(encode_png(f, 10), std::invalid_argument);
    if (detail::zlib_available())
        EXPECT_LT(encode_png(f, 9).size(), encode_png(f).size());
}

TEST(RenderSnapshot, AsciicastIsReproducible) {
    TextOptions t;
    t.mCols = 16;
    t.mRows = 8;
    SnapshotOptions s;
    s.mCastFrames = 4;
    s.mCastFps = 2.0;
    const std::string cast = encode_cast(cube_quads(), RenderOptions(), t, s);
    EXPECT_EQ(cast, encode_cast(cube_quads(), RenderOptions(), t, s));
    EXPECT_EQ(cast.rfind("{\"height\": 8, \"version\": 2, \"width\": 16}\n[0.000000, \"o\", \"", 0),
              0u);
    EXPECT_NE(cast.find("\n[1.500000, \"o\", \""), std::string::npos);
    EXPECT_EQ(std::count(cast.begin(), cast.end(), '\n'), 5);
}

TEST(RenderSnapshot, ChoosesTheFormByExtension) {
    EXPECT_THROW(write_snapshot("frame.vtu", cube_quads()), std::invalid_argument);
    TextOptions t;
    t.mEncoding = TextEncoding::Kitty;
    EXPECT_THROW(write_snapshot("frame.txt", cube_quads(), RenderOptions(), t),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Field rendering (v16.34.0): the pure helpers
// ---------------------------------------------------------------------------

TEST(RenderField, ScalesAreOddAndMonotone) {
    using detail::render_scale_forward;
    EXPECT_DOUBLE_EQ(render_scale_forward(RenderScale::Linear, 1.0, -3.5), -3.5);
    EXPECT_DOUBLE_EQ(render_scale_forward(RenderScale::Log, 1.0, 1000.0), 3.0);
    EXPECT_TRUE(std::isnan(render_scale_forward(RenderScale::Log, 1.0, 0.0)));
    EXPECT_TRUE(std::isnan(render_scale_forward(RenderScale::Log, 1.0, -1.0)));
    EXPECT_TRUE(std::isnan(render_scale_forward(RenderScale::Linear, 1.0, std::nan(""))));
    EXPECT_DOUBLE_EQ(render_scale_forward(RenderScale::Symlog, 1.0, 0.0), 0.0);
    EXPECT_DOUBLE_EQ(render_scale_forward(RenderScale::Symlog, 1.0, 9.0), 1.0);
    EXPECT_DOUBLE_EQ(render_scale_forward(RenderScale::Symlog, 1.0, -9.0), -1.0);
    double last = -1e300;
    for (double v = -100.0; v <= 100.0; v += 0.5) {
        const double t = render_scale_forward(RenderScale::Symlog, 0.5, v);
        EXPECT_GT(t, last);
        last = t;
    }
}

TEST(RenderField, PercentilesInterpolateBetweenRanks) {
    std::vector<double> v = {5, 1, 4, 2, 3};
    EXPECT_DOUBLE_EQ(detail::render_percentile(v, 0.0), 1.0);
    EXPECT_DOUBLE_EQ(detail::render_percentile(v, 100.0), 5.0);
    EXPECT_DOUBLE_EQ(detail::render_percentile(v, 50.0), 3.0);
    EXPECT_DOUBLE_EQ(detail::render_percentile(v, 25.0), 2.0);
    EXPECT_DOUBLE_EQ(detail::render_percentile(v, 10.0), 1.4);
    std::vector<double> none;
    EXPECT_TRUE(std::isnan(detail::render_percentile(none, 50.0)));
}

TEST(RenderField, LegendTicksAreRoundNumbers) {
    using detail::render_legend_ticks;
    EXPECT_EQ(render_legend_ticks(RenderScale::Linear, 1.0, 0.0, 1.0, 5),
              (std::vector<double>{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}));
    EXPECT_EQ(render_legend_ticks(RenderScale::Log, 1.0, 0.5, 2000.0, 5),
              (std::vector<double>{1.0, 10.0, 100.0, 1000.0}));
    const std::vector<double> sym = render_legend_ticks(RenderScale::Symlog, 1.0, -120.0, 50.0, 5);
    EXPECT_EQ(sym, (std::vector<double>{-100.0, -10.0, -1.0, 0.0, 1.0, 10.0}));
    EXPECT_TRUE(render_legend_ticks(RenderScale::Linear, 1.0, 2.0, 2.0, 5).empty());
}

TEST(RenderField, ContourOfATriangleIsOneSegmentAndConsistentAcrossAnEdge) {
    const double a[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    const double va[3] = {0.0, 1.0, 1.0};
    double seg[6];
    ASSERT_TRUE(detail::render_triangle_contour(a, va, 0.5, seg));
    EXPECT_DOUBLE_EQ(seg[0], 0.5);  // on edge 0-1
    EXPECT_DOUBLE_EQ(seg[3], 0.0);
    EXPECT_DOUBLE_EQ(seg[4], 0.5);  // on edge 2-0
    EXPECT_FALSE(detail::render_triangle_contour(a, va, 2.0, seg));
    // A level through a vertex: the corner counts as above, so both triangles of
    // a shared edge agree whether it is cut.
    const double vb[3] = {0.5, 1.0, 0.0};
    EXPECT_TRUE(detail::render_triangle_contour(a, vb, 0.5, seg));
    const double nan3[3] = {0.0, std::nan(""), 1.0};
    EXPECT_FALSE(detail::render_triangle_contour(a, nan3, 0.5, seg));
}

TEST(RenderField, PaletteHasTenDistinctColours) {
    std::set<std::uint32_t> seen;
    for (const auto& c : detail::render_category_palette())
        seen.insert((std::uint32_t(c[0]) << 16) | (std::uint32_t(c[1]) << 8) | c[2]);
    EXPECT_EQ(seen.size(), 10u);
}

// ---------------------------------------------------------------------------
// Field rendering: through render()
// ---------------------------------------------------------------------------

namespace {

// A unit square of n x n quads in z = 0 with a point array "x" equal to the
// x coordinate, "y" to y, and a cell array "id" of 0..n*n-1.
Mesh grid_square(int N) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i)
            pts.push_back({double(i) / N, double(j) / N, 0.0});
    std::vector<std::vector<std::int64_t>> quads;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const std::int64_t p = j * (N + 1) + i;
            quads.push_back({p, p + 1, p + N + 2, p + N + 1});
        }
    Mesh m = make_mesh(pts, "quad", quads);
    std::vector<std::vector<double>> xs, ys, ids;
    for (const auto& p : pts) {
        xs.push_back({p[0]});
        ys.push_back({p[1]});
    }
    for (int c = 0; c < N * N; ++c)
        ids.push_back({double(c)});
    m.AddPointData("x", mt::points_from(xs));
    m.AddPointData("y", mt::points_from(ys));
    m.AddCellData("id", {mt::points_from(ids)});
    return m;
}

RenderOptions top_view(int W = 64, int H = 64) {
    RenderOptions o;
    o.mWidth = W;
    o.mHeight = H;
    o.mView = "+z";
    o.mShading = RenderShading::None;
    return o;
}

std::size_t count_color(const Frame& rF, const RenderColor& rC) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < rF.mRgba.size(); i += 4)
        n += rF.mRgba[i] == rC[0] && rF.mRgba[i + 1] == rC[1] && rF.mRgba[i + 2] == rC[2] &&
                     rF.mRgba[i + 3] == 255
                 ? 1
                 : 0;
    return n;
}

NDArray cell_entries(const std::vector<std::int64_t>& rCells) {
    NDArray a(DType::Int64, {rCells.size()});
    for (std::size_t i = 0; i < rCells.size(); ++i)
        a.As<std::int64_t>()[i] = rCells[i];
    return a;
}

bool has_note(const Frame& rF, const std::string& rNeedle) {
    for (const std::string& n : rF.mNotes)
        if (n.find(rNeedle) != std::string::npos)
            return true;
    return false;
}

}  // namespace

TEST(RenderFieldFrame, PercentileClipEqualsTheExplicitRange) {
    // One huge outlier flattens the plot; clipping at 90 percent is the same
    // picture as giving the range that percentile names.
    Mesh m = grid_square(4);
    std::vector<std::vector<double>> v;
    for (int i = 0; i < 25; ++i)
        v.push_back({i == 24 ? 1000.0 : double(i)});
    m.AddPointData("u", mt::points_from(v));
    RenderOptions clip = top_view();
    clip.mColorBy = "u";
    clip.mClipHigh = 90.0;
    const Frame a = render(m, clip);
    RenderOptions explicit_range = top_view();
    explicit_range.mColorBy = "u";
    explicit_range.mVMin = a.mVMin;
    explicit_range.mVMax = a.mVMax;
    const Frame b = render(m, explicit_range);
    EXPECT_LT(a.mVMax, 100.0);
    EXPECT_EQ(a.mRgba, b.mRgba);
    clip.mClipLow = 90.0;
    clip.mClipHigh = 10.0;
    EXPECT_THROW(render(m, clip), std::invalid_argument);
}

TEST(RenderFieldFrame, SymmetricRangePutsZeroAtTheMidpoint) {
    Mesh m = grid_square(4);
    std::vector<std::vector<double>> v;
    for (int i = 0; i < 25; ++i)
        v.push_back({double(i) - 4.0});  // -4 .. 20
    m.AddPointData("u", mt::points_from(v));
    RenderOptions o = top_view();
    o.mColorBy = "u";
    o.mSymmetric = true;
    o.mCmap = "coolwarm";
    const Frame f = render(m, o);
    // A face's value is the mean of its corners, so the extremes are not -4 and
    // 20; what matters is the range is centred on zero.
    EXPECT_GT(f.mVMax, 10.0);
    EXPECT_DOUBLE_EQ(f.mVMin, -f.mVMax);
}

TEST(RenderFieldFrame, LogScaleRefusesANonPositiveRangeAndReportsIt) {
    Mesh m = grid_square(2);  // four cells: cell data are drawn exactly
    m.AddCellData("c", {mt::points_from({{1}, {10}, {100}, {1000}})});
    RenderOptions o = top_view();
    o.mColorBy = "c";
    o.mScale = RenderScale::Log;
    o.mColorbar = true;
    const Frame f = render(m, o);
    EXPECT_DOUBLE_EQ(f.mVMin, 1.0);
    EXPECT_DOUBLE_EQ(f.mVMax, 1000.0);
    EXPECT_TRUE(has_note(f, "log"));
    EXPECT_TRUE(has_note(f, "ticks: 1, 10, 100, 1000"));
    o.mVMin = 0.0;
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(RenderFieldFrame, NonFiniteValuesAreCountedInTheNotes) {
    Mesh m = grid_square(2);
    m.AddPointData("u", mt::points_from({{1}, {2}, {3}, {std::nan("")}, {5}, {6}, {7}, {8}, {9}}));
    RenderOptions o = top_view();
    o.mColorBy = "u";
    const Frame f = render(m, o);
    EXPECT_TRUE(has_note(f, "no finite value"));
}

TEST(RenderFieldFrame, IsolinesOfXAreEquallySpacedVerticalLines) {
    // 8 x 8 quads in the unit square, seen from +z, 80 pixels wide: the points
    // x = 1/8 ... 7/8 are fixed by the grid, so the contour levels 1/4 .. 3/4
    // (three of them) fall on columns a quarter, a half and three quarters of the way
    // across the drawn square.
    Mesh m = grid_square(8);
    RenderOptions o = top_view(80, 80);
    o.mColorBy = "x";
    o.mIsolines = 3;
    o.mIsoColor = {255, 0, 255, 255};
    o.mCmap = "grey";
    const Frame f = render(m, o);
    std::vector<int> columns;
    for (int x = 0; x < f.mWidth; ++x) {
        int hits = 0;
        for (int y = 0; y < f.mHeight; ++y) {
            const std::uint8_t* p = f.mRgba.data() + (std::size_t(y) * f.mWidth + x) * 4;
            hits += p[0] == 255 && p[1] == 0 && p[2] == 255 ? 1 : 0;
        }
        if (hits > f.mHeight / 2)
            columns.push_back(x);
    }
    // Each line is one or two pixels wide; collapse neighbours.
    std::vector<int> centres;
    for (std::size_t i = 0; i < columns.size(); ++i)
        if (i == 0 || columns[i] != columns[i - 1] + 1)
            centres.push_back(columns[i]);
    ASSERT_EQ(centres.size(), 3u);
    EXPECT_NEAR(centres[1] - centres[0], centres[2] - centres[1], 1.5);
    // They are vertical: no row of the lines is longer than a few pixels.
    o.mIsoLevels = {0.5};
    o.mIsolines = 0;
    EXPECT_GT(count_color(render(m, o), {255, 0, 255, 255}), 20u);
    // Contours of cell data are refused by name.
    RenderOptions bad = top_view();
    bad.mColorBy = "id";
    bad.mIsolines = 2;
    EXPECT_THROW(render(m, bad), std::invalid_argument);
}

TEST(RenderFieldFrame, ColoursByRegionsWithAKeyThatListsEmptyOnesToo) {
    Mesh m = grid_square(2);  // four cells
    m.AddRegion(Region("left", RegionKind::Cell, cell_entries({0, 2})));
    m.AddRegion(Region("right", RegionKind::Cell, cell_entries({1, 3})));
    m.AddRegion(Region("none", RegionKind::Cell, cell_entries({})));
    RenderOptions o = top_view();
    o.mColorRegions = true;
    const Frame f = render(m, o);
    EXPECT_TRUE(has_note(f, "region left: #0072b2 (2 cells)"));
    // Regions come back in (kind, name) order, which fixes the palette order.
    EXPECT_TRUE(has_note(f, "region none: #e69f00 (0 cells)"));
    EXPECT_TRUE(has_note(f, "region right: #009e73 (2 cells)"));
    EXPECT_GT(count_color(f, {0, 114, 178, 255}), 0u);
    EXPECT_GT(count_color(f, {0, 158, 115, 255}), 0u);
    EXPECT_EQ(count_color(f, {230, 159, 0, 255}), 0u);  // the empty region paints nothing
    o.mCategoryEdges = true;
    o.mEdgeColor = {255, 255, 255, 255};
    EXPECT_GT(count_color(render(m, o), {255, 255, 255, 255}), 0u);
    RenderOptions none = top_view();
    none.mColorRegions = true;
    EXPECT_THROW(render(grid_square(2), none), std::invalid_argument);
}

TEST(RenderFieldFrame, CategoricalDataGetsAKeyWithCounts) {
    Mesh m = grid_square(2);
    m.AddCellData("mat", {mt::points_from({{7}, {7}, {9}, {9}})});
    RenderOptions o = top_view();
    o.mColorBy = "mat";
    o.mCategorical = true;
    const Frame f = render(m, o);
    EXPECT_TRUE(has_note(f, "mat = 7: #0072b2 (2 drawn faces)"));
    EXPECT_TRUE(has_note(f, "mat = 9: #e69f00 (2 drawn faces)"));
    EXPECT_FALSE(f.mColored);
}

TEST(RenderFieldFrame, ExpressionsAndTensorReductionsColourLikeArrays) {
    Mesh m = grid_square(2);
    m.AddPointData("u", mt::points_from({{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}}));
    RenderOptions o = top_view();
    o.mExpr = "u * 2";
    const Frame f = render(m, o);
    // Face values are corner means: the first quad's corners 0, 1, 3, 4 give 2,
    // the last's 4, 5, 7, 8 give 6; doubled by the expression.
    EXPECT_DOUBLE_EQ(f.mVMin, 4.0);
    EXPECT_DOUBLE_EQ(f.mVMax, 12.0);
    EXPECT_TRUE(has_note(f, "u * 2: 4 .. 12"));
    o.mColorBy = "u";  // both: refused
    EXPECT_THROW(render(m, o), std::invalid_argument);
    o.mExpr = "nope + 1";
    o.mColorBy.clear();
    EXPECT_THROW(render(m, o), std::invalid_argument);

    // A uniaxial tension of 5 everywhere: von Mises 5, hydrostatic 5/3, and
    // principal values 0, 0, 5.
    std::vector<std::vector<double>> tensor(9, {5, 0, 0, 0, 0, 0});
    m.AddPointData("stress", mt::points_from(tensor));
    RenderOptions r = top_view();
    r.mColorBy = "stress";
    r.mReduce = "mises";
    EXPECT_NEAR(render(m, r).mVMax, 5.0, 1e-12);
    r.mReduce = "hydrostatic";
    EXPECT_NEAR(render(m, r).mVMax, 5.0 / 3.0, 1e-12);
    r.mReduce = "principal";
    EXPECT_NEAR(render(m, r).mVMax, 5.0, 1e-12);
    r.mComponent = 0;
    EXPECT_NEAR(render(m, r).mVMax, 0.0, 1e-12);
    r.mReduce = "nonsense";
    EXPECT_THROW(render(m, r), std::invalid_argument);
}

TEST(RenderFieldFrame, WarpMovesThePointsAndDrawsTheOutline) {
    Mesh m = grid_square(2);
    std::vector<std::vector<double>> d;
    for (int i = 0; i < 9; ++i)
        d.push_back({0.0, 0.0, 0.5});
    m.AddPointData("u", mt::points_from(d));
    RenderOptions base = top_view();
    base.mView = "+x";  // from the side a lift in z is visible
    const Frame before = render(m, base);
    RenderOptions o = base;
    o.mWarp = "u";
    o.mWarpOutline = true;
    o.mOutlineColor = {255, 0, 0, 255};
    const Frame after = render(m, o);
    EXPECT_NE(before.mRgba, after.mRgba);
    EXPECT_GT(count_color(after, {255, 0, 0, 255}), 0u);
    EXPECT_TRUE(has_note(after, "warp: u x 1, undeformed outline"));
    o.mWarp = "nope";
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(RenderFieldFrame, VectorsDrawArrowsInTheirColour) {
    Mesh m = grid_square(4);
    std::vector<std::vector<double>> v(25, {0.0, 1.0, 0.0});
    m.AddPointData("v", mt::points_from(v));
    RenderOptions o = top_view(120, 120);
    o.mVectors = "v";
    o.mVectorCount = 9;
    o.mVectorColor = {255, 0, 0, 255};
    const Frame f = render(m, o);
    EXPECT_GT(count_color(f, {255, 0, 0, 255}), 30u);
    EXPECT_TRUE(has_note(f, "vectors: v, 9 arrows"));
    o.mVectorLength = 0.05;
    EXPECT_GT(count_color(render(m, o), {255, 0, 0, 255}), 10u);
}

TEST(RenderFieldFrame, ADiagnosticShowsAFlippedTriangle) {
    // Two triangles on a square; the second is wound the other way, so seen
    // from +z one is a front face and one a back face.
    Mesh m =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}, {0, 2, 3}});
    RenderOptions o = top_view();
    o.mDiagnostic = RenderDiagnostic::Orientation;
    Frame f = render(m, o);
    EXPECT_GT(count_color(f, {70, 130, 230, 255}), 0u);
    EXPECT_EQ(count_color(f, {230, 140, 40, 255}), 0u);
    Mesh flipped =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}, {0, 3, 2}});
    f = render(flipped, o);
    EXPECT_GT(count_color(f, {70, 130, 230, 255}), 0u);
    EXPECT_GT(count_color(f, {230, 140, 40, 255}), 0u);
    EXPECT_TRUE(has_note(f, "orientation:"));
}

TEST(RenderFieldFrame, OtherDiagnostics) {
    Mesh m =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle", {{0, 1, 2}, {0, 2, 3}});
    RenderOptions o = top_view();
    o.mDiagnostic = RenderDiagnostic::FreeEdges;
    Frame f = render(m, o);
    EXPECT_GT(count_color(f, {235, 140, 0, 255}), 0u);  // the open boundary
    EXPECT_TRUE(has_note(f, "free edges: 4 open"));

    o = top_view();
    o.mDiagnostic = RenderDiagnostic::EdgeLength;
    f = render(m, o);
    EXPECT_TRUE(f.mColored);
    EXPECT_TRUE(has_note(f, "edge length:"));

    o = top_view();
    o.mDiagnostic = RenderDiagnostic::Quality;
    o.mQualityMetric = "scaled_jacobian";
    f = render(m, o);
    EXPECT_TRUE(f.mColored);
    o.mQualityMetric = "bogus";
    EXPECT_THROW(render(m, o), std::invalid_argument);
    o.mQualityMetric.clear();
    EXPECT_THROW(render(m, o), std::invalid_argument);

    // A triangle wound against its neighbour's orientation is "inverted" only
    // by signed area in a 3-D sense, so just check the flag path runs and keys.
    o = top_view();
    o.mDiagnostic = RenderDiagnostic::Degenerate;
    f = render(m, o);
    EXPECT_TRUE(has_note(f, "degenerate: 0 drawn faces in red"));
    o.mColorBy = "x";
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(RenderFieldFrame, ColoursOfTheOriginalFieldAreUnchangedByTheNewOptions) {
    // The v16.33.0 behaviour is the default of v16.34.0: same bytes.
    Mesh m = grid_square(4);
    RenderOptions a = top_view();
    a.mColorBy = "x";
    a.mColorbar = true;
    RenderOptions b = a;
    b.mScale = RenderScale::Linear;
    b.mSymmetric = false;
    EXPECT_EQ(render(m, a).mRgba, render(m, b).mRgba);
}

// ---------------------------------------------------------------------------
// Cut-aways (v16.36.0, ABI 25): planes that clip the drawn geometry
// ---------------------------------------------------------------------------

namespace {

RenderCutaway keep_beyond(int Axis, double At, double Sign = 1.0) {
    RenderCutaway c;
    c.mPoint = {0.0, 0.0, 0.0};
    c.mPoint[static_cast<std::size_t>(Axis)] = At;
    c.mNormal = {0.0, 0.0, 0.0};
    c.mNormal[static_cast<std::size_t>(Axis)] = Sign;
    return c;
}

RenderOptions side_view(const char* pView, int Size = 64) {
    RenderOptions o;
    o.mWidth = Size;
    o.mHeight = Size;
    o.mView = pView;
    o.mShading = RenderShading::None;
    o.mBackground = {0, 0, 0, 255};
    return o;
}

std::int64_t centre_id(const Frame& rFrame) {
    return rFrame.mCellIds[static_cast<std::size_t>(rFrame.mHeight / 2 * rFrame.mWidth +
                                                    rFrame.mWidth / 2)];
}

RenderColor centre_rgba(const Frame& rFrame) {
    const std::size_t at = static_cast<std::size_t>(rFrame.mHeight / 2 * rFrame.mWidth +
                                                    rFrame.mWidth / 2) *
                           4;
    return {rFrame.mRgba[at], rFrame.mRgba[at + 1], rFrame.mRgba[at + 2], rFrame.mRgba[at + 3]};
}

std::size_t covered(const Frame& rFrame, std::int64_t Id) {
    return static_cast<std::size_t>(
        std::count(rFrame.mCellIds.begin(), rFrame.mCellIds.end(), Id));
}

}  // namespace

TEST(RenderCutaway, ACutThroughACubeShowsTheTintedInsideOfTheFarFace) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("-x");  // the camera sits on the -x side
    EXPECT_EQ(centre_id(render(m, o)), 4);  // the x = 0 face, uncut
    o.mCutaways = {keep_beyond(0, 0.5)};
    o.mCutawayTint = {10, 200, 30, 255};
    const Frame cut = render(m, o);
    EXPECT_EQ(centre_id(cut), 5);  // now the x = 1 face, seen from inside
    EXPECT_EQ(centre_rgba(cut), (RenderColor{10, 200, 30, 255}));
    EXPECT_EQ(covered(cut, 4), 0u);  // the near face is gone
    EXPECT_TRUE(std::any_of(cut.mNotes.begin(), cut.mNotes.end(), [](const std::string& s) {
        return s.rfind("cutaway: 1 plane", 0) == 0;
    }));
}

TEST(RenderCutaway, FrontFacesAreNotTinted) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("+x");  // the camera sits on the +x side
    o.mCutaways = {keep_beyond(0, 0.5)};
    const Frame cut = render(m, o);
    EXPECT_EQ(centre_id(cut), 5);  // the x = 1 face, seen from outside
    EXPECT_EQ(centre_rgba(cut), (RenderColor{200, 197, 189, 255}));
}

TEST(RenderCutaway, TheKeptSideIsTheOneTheNormalPointsTo) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("+x");
    o.mCutaways = {keep_beyond(0, 0.5, -1.0)};  // keeps x <= 0.5: the x = 1 face is gone
    const Frame cut = render(m, o);
    EXPECT_EQ(covered(cut, 5), 0u);
    EXPECT_EQ(centre_id(cut), 4);  // the x = 0 face, from inside, tinted
}

TEST(RenderCutaway, TwoPlanesLeaveAQuarter) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("+z");
    const std::size_t whole = covered(render(m, o), 1);  // the z = 1 face
    ASSERT_GT(whole, 0u);
    o.mCutaways = {keep_beyond(0, 0.5), keep_beyond(1, 0.5)};
    const std::size_t quarter = covered(render(m, o), 1);
    EXPECT_GT(quarter, whole / 5);
    EXPECT_LT(quarter, whole * 3 / 10);
}

TEST(RenderCutaway, AFaceCrossingThePlaneKeepsItsIdAndColour) {
    const Mesh m = grid_square(4);
    RenderOptions o = top_view();
    o.mColorBy = "x";
    const Frame whole = render(m, o);
    o.mCutaways = {keep_beyond(0, 0.30)};
    const Frame cut = render(m, o);
    // Cell values are corner means, so a cell cut in two keeps one colour; no
    // pixel appears where there was none.
    for (std::size_t i = 0; i < cut.mCellIds.size(); ++i)
        if (cut.mCellIds[i] >= 0) {
            EXPECT_EQ(cut.mCellIds[i], whole.mCellIds[i]);
        }
    EXPECT_LT(covered(cut, -1), cut.mCellIds.size());
    EXPECT_GT(covered(cut, -1), covered(whole, -1));
}

TEST(RenderCutaway, LinesAndPointsAreClippedToo) {
    const Mesh lines = make_mesh({{0, 0, 0}, {1, 0, 0}}, "line", {{0, 1}});
    RenderOptions o = top_view();
    o.mLineColor = {255, 0, 0, 255};
    const std::size_t whole = covered(render(lines, o), 0);
    ASSERT_GT(whole, 4u);
    o.mCutaways = {keep_beyond(0, 0.5)};
    const std::size_t half = covered(render(lines, o), 0);
    EXPECT_GT(half, whole / 3);
    EXPECT_LT(half, whole * 2 / 3);
    o.mCutaways = {keep_beyond(0, 2.0)};
    EXPECT_EQ(covered(render(lines, o), 0), 0u);  // nothing on the kept side

    const Mesh cloud = make_mesh({{0, 0, 0}, {1, 0, 0}, {0.2, 0.5, 0}}, "vertex", {{0}, {1}, {2}});
    RenderOptions p = top_view(64, 64);
    p.mPointRadius = 2.0;
    p.mCutaways = {keep_beyond(0, 0.5)};
    const Frame points = render(cloud, p);
    EXPECT_GT(covered(points, 1), 0u);
    EXPECT_EQ(covered(points, 0), 0u);
    EXPECT_EQ(covered(points, 2), 0u);
}

TEST(RenderCutaway, ASceneDrawnWithACutEqualsRenderWithTheCut) {
    const Mesh m = cube_quads();
    RenderOptions base = side_view("-x");
    base.mShading = RenderShading::Smooth;
    base.mEdges = RenderEdges::All;
    base.mBackground = {0, 0, 0, 0};
    const RenderScene scene = prepare_render(m, base);
    RenderOptions cut = base;
    cut.mCutaways = {keep_beyond(0, 0.35), keep_beyond(2, 0.25)};
    cut.mAzimuth = 20.0;
    cut.mView.clear();
    const Frame direct = render(m, cut);
    const Frame drawn = render_scene(scene, cut);
    EXPECT_EQ(direct.mRgba, drawn.mRgba);
    EXPECT_EQ(direct.mCellIds, drawn.mCellIds);
    // And the cut is a draw-time choice: the same scene draws uncut again.
    EXPECT_EQ(render(m, base).mRgba, render_scene(scene, base).mRgba);
    EXPECT_NE(direct.mRgba, render_scene(scene, base).mRgba);
}

TEST(RenderCutaway, ThePerspectiveCameraAndFieldsWork) {
    const Mesh m = grid_square(4);
    RenderOptions o = top_view();
    o.mProjection = RenderProjection::Perspective;
    o.mColorBy = "x";
    o.mIsolines = 3;
    o.mVectors = "x";
    o.mCutaways = {keep_beyond(1, 0.4)};
    EXPECT_NO_THROW(render(m, o));
}

TEST(RenderCutaway, ClippingEverythingIsABlankFrameNotAnError) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("+z");
    o.mCutaways = {keep_beyond(2, 10.0)};
    const Frame f = render(m, o);
    EXPECT_EQ(covered(f, -1), f.mCellIds.size());
}

TEST(RenderCutaway, BadPlanesAreRefused) {
    const Mesh m = cube_quads();
    RenderOptions o = side_view("+z");
    o.mCutaways = {keep_beyond(0, 0.1), keep_beyond(1, 0.1), keep_beyond(2, 0.1)};
    EXPECT_THROW(render(m, o), std::invalid_argument);
    RenderCutaway flat;
    flat.mNormal = {0.0, 0.0, 0.0};
    o.mCutaways = {flat};
    EXPECT_THROW(render(m, o), std::invalid_argument);
    RenderCutaway nan_plane;
    nan_plane.mPoint = {std::nan(""), 0.0, 0.0};
    o.mCutaways = {nan_plane};
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(RenderCutaway, WithNoPlanesTheBytesAreUnchanged) {
    const Mesh m = cube_quads();
    RenderOptions a = side_view("+x");
    a.mShading = RenderShading::Smooth;
    RenderOptions b = a;
    b.mCutaways.clear();
    b.mCutawayTint = {1, 2, 3, 255};  // the tint matters only when something is cut
    EXPECT_EQ(render(m, a).mRgba, render(m, b).mRgba);
}

// ---------------------------------------------------------------------------
// Streamlines
// ---------------------------------------------------------------------------

namespace {

NDArray swirl_field(const Mesh& rMesh) {
    std::vector<std::vector<double>> v;
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
        const double x = detail::read_double(pts, i * dim + 0);
        const double y = detail::read_double(pts, i * dim + 1);
        v.push_back({-(y - 0.5), x - 0.5, 0.0});
    }
    return mt::points_from(v);
}

std::vector<double> flat_field(const NDArray& rField) {
    std::vector<double> out(rField.Size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = detail::read_double(rField, i);
    return out;
}

std::vector<double> flat_points(const Mesh& rMesh) {
    std::vector<double> out;
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            out.push_back(k < dim ? detail::read_double(pts, i * dim + k) : 0.0);
    return out;
}

}  // namespace

TEST(Streamlines, ASwirlOverAGridFollowsCircles) {
    Mesh m = grid_square(8);
    const std::vector<double> xyz = flat_points(m);
    const std::vector<double> vec = flat_field(swirl_field(m));
    detail::StreamlineOptions o;
    o.mSeeds = 12;
    o.mLength = 0.4;
    const detail::Streamlines lines = detail::trace_streamlines(m, xyz.data(), vec.data(), o);
    EXPECT_EQ(lines.mDim, 2);
    EXPECT_GE(lines.NumLines(), 8u);
    for (std::size_t l = 0; l < lines.NumLines(); ++l) {
        const std::size_t first = lines.mStart[l];
        auto radius = [&](std::size_t k) {
            return std::hypot(lines.mXyz[3 * k] - 0.5, lines.mXyz[3 * k + 1] - 0.5);
        };
        for (std::size_t k = first; k < lines.mStart[l + 1]; ++k)
            EXPECT_NEAR(radius(k), radius(first), 2e-3) << "line " << l;
    }
}

TEST(Streamlines, TheSameTraceComesBackEveryTime) {
    Mesh m = grid_square(6);
    const std::vector<double> xyz = flat_points(m);
    const std::vector<double> vec = flat_field(swirl_field(m));
    const detail::Streamlines a =
        detail::trace_streamlines(m, xyz.data(), vec.data(), detail::StreamlineOptions{});
    const detail::Streamlines b =
        detail::trace_streamlines(m, xyz.data(), vec.data(), detail::StreamlineOptions{});
    EXPECT_EQ(a.mXyz, b.mXyz);
    EXPECT_EQ(a.mStart, b.mStart);
}

TEST(Streamlines, AUniformFieldCrossesAHexahedronToItsFaces) {
    Mesh m = mt::hex_mesh();
    const std::vector<double> xyz = flat_points(m);
    std::vector<double> vec;
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        vec.insert(vec.end(), {1.0, 0.0, 0.0});
    detail::StreamlineOptions o;
    o.mLength = 2.0;
    const detail::Streamlines lines = detail::trace_streamlines(m, xyz.data(), vec.data(), o);
    EXPECT_EQ(lines.mDim, 3);
    ASSERT_GT(lines.NumLines(), 0u);
    for (std::size_t l = 0; l < lines.NumLines(); ++l) {
        double lo = 1e9;
        double hi = -1e9;
        for (std::size_t k = lines.mStart[l]; k < lines.mStart[l + 1]; ++k) {
            EXPECT_NEAR(lines.mXyz[3 * k + 1], lines.mXyz[3 * lines.mStart[l] + 1], 1e-9);
            EXPECT_NEAR(lines.mXyz[3 * k + 2], lines.mXyz[3 * lines.mStart[l] + 2], 1e-9);
            lo = std::min(lo, lines.mXyz[3 * k]);
            hi = std::max(hi, lines.mXyz[3 * k]);
        }
        EXPECT_LT(lo, 0.02);
        EXPECT_GT(hi, 0.98);
    }
}

TEST(Streamlines, ATiltedSurfaceKeepsTheLinesOnIt) {
    // The plane z = x, with a field along it.
    Mesh m = make_mesh({{0, 0, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 0}}, "triangle",
                       {{0, 1, 2}, {0, 2, 3}});
    const std::vector<double> xyz = flat_points(m);
    std::vector<double> vec;
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        vec.insert(vec.end(), {1.0, 0.0, 1.0});
    detail::StreamlineOptions o;
    o.mSeeds = 4;
    o.mLength = 1.0;
    const detail::Streamlines lines = detail::trace_streamlines(m, xyz.data(), vec.data(), o);
    EXPECT_EQ(lines.mDim, 2);
    ASSERT_GT(lines.NumLines(), 0u);
    for (std::size_t k = 0; k < lines.mXyz.size() / 3; ++k)
        EXPECT_NEAR(lines.mXyz[3 * k + 2], lines.mXyz[3 * k], 1e-9);
}

TEST(Streamlines, ACellTypeWithNothingToFollowIsNamed) {
    Mesh m = mt::line_mesh();
    const std::vector<double> xyz = flat_points(m);
    const std::vector<double> vec(3 * m.NumPoints(), 1.0);
    try {
        detail::trace_streamlines(m, xyz.data(), vec.data(), detail::StreamlineOptions{});
        FAIL() << "a mesh of lines has no streamlines";
    } catch (const std::invalid_argument& rErr) {
        EXPECT_NE(std::string(rErr.what()).find("line"), std::string::npos);
        EXPECT_NE(std::string(rErr.what()).find("streamlines need"), std::string::npos);
    }
}

TEST(Streamlines, TheRendererDrawsThemInTheirColour) {
    Mesh m = grid_square(8);
    m.AddPointData("v", swirl_field(m));
    RenderOptions o = top_view(120, 120);
    EXPECT_EQ(count_color(render(m, o), {240, 80, 160, 255}), 0u);
    o.mStreamlines = "v";
    o.mStreamColor = {255, 0, 0, 255};
    const Frame f = render(m, o);
    EXPECT_GT(count_color(f, {255, 0, 0, 255}), 100u);
    EXPECT_TRUE(has_note(f, "streamlines: v,"));
    o.mStreamSeeds = 5;
    const Frame fewer = render(m, o);
    EXPECT_LT(count_color(fewer, {255, 0, 0, 255}), count_color(f, {255, 0, 0, 255}));
    o.mStreamlines = "nope";
    EXPECT_THROW(render(m, o), std::invalid_argument);
    o.mStreamlines = "v";
    o.mStreamSeeds = 0;
    EXPECT_THROW(render(m, o), std::invalid_argument);
    o.mStreamSeeds = 5;
    o.mStreamLength = 0.0;
    EXPECT_THROW(render(m, o), std::invalid_argument);
}

TEST(Streamlines, OffByDefaultTheBytesAreUnchanged) {
    Mesh m = grid_square(4);
    m.AddPointData("v", swirl_field(m));
    RenderOptions a = top_view(80, 80);
    RenderOptions b = a;
    b.mStreamColor = {1, 2, 3, 255};
    b.mStreamSeeds = 7;  // nothing is drawn without a named array
    EXPECT_EQ(render(m, a).mRgba, render(m, b).mRgba);
}

// ---------------------------------------------------------------------------
// The synthwave theme
// ---------------------------------------------------------------------------

namespace {

RenderColor pixel(const Frame& rF, int X, int Y) {
    const std::size_t at = (static_cast<std::size_t>(Y) * rF.mWidth + X) * 4;
    return {rF.mRgba[at], rF.mRgba[at + 1], rF.mRgba[at + 2], rF.mRgba[at + 3]};
}

RenderOptions themed(int W = 64, int H = 64) {
    RenderOptions o = top_view(W, H);
    o.mTheme = RenderTheme::Synthwave;
    return o;
}

}  // namespace

TEST(RenderTheme, ABandedSunsetFillsTheTransparentBackgroundOnly) {
    const Mesh m = grid_square(2);
    RenderOptions o = themed(64, 96);
    o.mZoom = 0.4;  // leave the corners empty
    const Frame f = render(m, o);
    EXPECT_EQ(pixel(f, 0, 0), (RenderColor{20, 8, 60, 255}));  // the top of the sky
    EXPECT_EQ(pixel(f, 63, 0), pixel(f, 0, 0));
    EXPECT_EQ(pixel(f, 0, 95)[3], 255);
    // Banded: eight bands of sky and four of floor, not a smooth ramp.
    std::set<std::array<std::uint8_t, 3>> colours;
    for (int y = 0; y < 96; ++y) {
        const RenderColor c = pixel(f, 0, y);
        colours.insert({c[0], c[1], c[2]});
    }
    EXPECT_LE(colours.size(), 12u);
    EXPECT_GE(colours.size(), 8u);
    // The warm end sits at the horizon, above the floor.
    EXPECT_GT(pixel(f, 0, 55)[0], pixel(f, 0, 5)[0]);
    EXPECT_GT(pixel(f, 0, 55)[1], pixel(f, 0, 5)[1]);
    // An opaque background of your own wins.
    o.mBackground = {1, 2, 3, 255};
    EXPECT_EQ(pixel(render(m, o), 0, 0), (RenderColor{1, 2, 3, 255}));
}

TEST(RenderTheme, TheThemeColoursTheFacesAndEdgesUnlessYouSetThem) {
    const Mesh m = grid_square(2);
    RenderOptions o = themed();
    const Frame f = render(m, o);
    EXPECT_EQ(pixel(f, 32, 32), (RenderColor{96, 56, 190, 255}));  // neon violet, unlit
    o.mFillColor = {10, 20, 30, 255};
    EXPECT_EQ(pixel(render(m, o), 32, 32), (RenderColor{10, 20, 30, 255}));
    o = themed();
    o.mEdges = RenderEdges::All;
    EXPECT_GT(count_color(render(m, o), {0, 230, 255, 255}), 20u);  // cyan edges
}

TEST(RenderTheme, TheGridFloorNeedsATheme) {
    const Mesh m = grid_square(2);
    RenderOptions o = top_view();
    o.mGridFloor = true;
    EXPECT_THROW(render(m, o), std::invalid_argument);
    o.mTheme = RenderTheme::Synthwave;
    EXPECT_NO_THROW(render(m, o));
}

TEST(RenderTheme, TheGridFloorStepsWithTheBeat) {
    const Mesh m = grid_square(2);
    RenderOptions o = themed(96, 96);
    o.mZoom = 0.3;
    const Frame plain = render(m, o);
    o.mGridFloor = true;
    const Frame grid = render(m, o);
    EXPECT_NE(plain.mRgba, grid.mRgba);
    EXPECT_GT(count_color(grid, {255, 50, 200, 255}), 50u);  // magenta on an even beat
    o.mThemePhase = 1;
    const Frame next = render(m, o);
    EXPECT_GT(count_color(next, {0, 220, 255, 255}), 50u);  // cyan on an odd one
    EXPECT_EQ(count_color(next, {255, 50, 200, 255}), 0u);
    o.mThemePhase = 8;  // the scroll comes round after eight beats
    EXPECT_EQ(render(m, o).mRgba, grid.mRgba);
}

TEST(RenderTheme, ThePostProcessesChangeTheFrameAndAreDeterministic) {
    const Mesh m = grid_square(4);
    RenderOptions base = themed(120, 80);
    base.mColorBy = "x";
    base.mCmap = "synthwave";
    const Frame plain = render(m, base);
    for (int which = 0; which < 3; ++which) {
        RenderOptions o = base;
        (which == 0 ? o.mBloom : which == 1 ? o.mFringe : o.mScanlines) = true;
        const Frame a = render(m, o);
        EXPECT_NE(a.mRgba, plain.mRgba) << which;
        EXPECT_EQ(a.mRgba, render(m, o).mRgba) << which;
        EXPECT_EQ(a.mCellIds, plain.mCellIds) << which;  // pixels change, picking does not
    }
    // Scanlines dim every other row and leave the rest.
    RenderOptions o = base;
    o.mScanlines = true;
    const Frame lines = render(m, o);
    for (int x : {10, 60, 100}) {
        EXPECT_EQ(pixel(lines, x, 20), pixel(plain, x, 20));
        const RenderColor a = pixel(plain, x, 21);
        const RenderColor b = pixel(lines, x, 21);
        EXPECT_EQ(b[0], a[0] - a[0] / 4);
        EXPECT_EQ(b[3], 255);
    }
}

TEST(RenderTheme, WithEverythingOffTheBytesAreUnchanged) {
    const Mesh m = grid_square(4);
    RenderOptions a = top_view(80, 60);
    a.mColorBy = "x";
    RenderOptions b = a;
    b.mThemePhase = 5;  // the phase matters only with a theme
    EXPECT_EQ(render(m, a).mRgba, render(m, b).mRgba);
}
