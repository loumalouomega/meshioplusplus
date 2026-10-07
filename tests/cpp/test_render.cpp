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
#include <cstdint>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/operations/render.hpp"
#include "mesh_fixtures.hpp"

// Project includes (private, not installed)
#include "../../src/cpp/src/detail/png_write.hpp"
#include "../../src/cpp/src/detail/raster.hpp"

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
