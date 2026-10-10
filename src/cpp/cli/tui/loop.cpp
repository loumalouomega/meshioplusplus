// SPDX-License-Identifier: MIT
/// @file loop.cpp
/// @brief Implementation of loop.hpp.

#include "loop.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <set>
#include <utility>

#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "session_file.hpp"

namespace meshioplusplus::cli::tui {

namespace {

constexpr double kZoomStep = 1.1;
constexpr double kPanStep = 0.05;
constexpr double kMinZoom = 0.05;
constexpr double kMaxZoom = 200.0;
constexpr int kMaxProbeRows = 6;

// The named views, the same camera directions `RenderOptions::mView` selects
// (the library keeps its own copy private); the loop turns a drag into an
// angle change from wherever a preset left the camera.
struct Preset {
    const char* mName;
    double mAzimuth;
    double mElevation;
};
const Preset kPresets[] = {{"iso", 45.0, 35.264389682754654}, {"+x", 0.0, 0.0},
                           {"-x", 180.0, 0.0},                {"+y", 90.0, 0.0},
                           {"-y", -90.0, 0.0},                {"+z", -90.0, 90.0},
                           {"-z", -90.0, -90.0}};

// A named view as angles on `rOptions`, clearing the name.
void apply_view_name(RenderOptions& rOptions) {
    if (rOptions.mView.empty())
        return;
    for (const Preset& preset : kPresets)
        if (rOptions.mView == preset.mName) {
            rOptions.mAzimuth = preset.mAzimuth;
            rOptions.mElevation = preset.mElevation;
            rOptions.mView.clear();
            return;
        }
    throw std::invalid_argument("unknown view '" + rOptions.mView +
                                "' (expected iso, +x, -x, +y, -y, +z or -z)");
}

bool is_cell_encoding(TextEncoding Encoding) {
    return Encoding != TextEncoding::Kitty && Encoding != TextEncoding::ITerm2 &&
           Encoding != TextEncoding::Sixel;
}

// Pixels per terminal cell, horizontally and vertically.
std::array<int, 2> cell_pixels(TextEncoding Encoding) {
    switch (Encoding) {
        case TextEncoding::Quadrant:
            return {2, 2};
        case TextEncoding::Sextant:
            return {2, 3};
        case TextEncoding::Braille:
            return {2, 4};
        case TextEncoding::Ascii:
            return {1, 1};
        default:
            return {1, 2};
    }
}

// `Text` cut or padded to exactly `Cols` columns (one column per code point:
// the status and notes are ASCII, and a wide glyph in a title only shifts the
// padding).
std::string fit_columns(const std::string& rText, int Cols) {
    std::string out;
    int used = 0;
    for (std::size_t i = 0; i < rText.size() && used < Cols;) {
        const unsigned char c = static_cast<unsigned char>(rText[i]);
        const std::size_t len = c < 0x80 ? 1 : c >= 0xF0 ? 4 : c >= 0xE0 ? 3 : 2;
        out.append(rText, i, std::min(len, rText.size() - i));
        i += len;
        ++used;
    }
    if (used < Cols)
        out.append(static_cast<std::size_t>(Cols - used), ' ');
    return out;
}

std::string at_row(int Row) {
    return "\x1b[" + std::to_string(Row) + ";1H";
}

std::string number(double Value, int Decimals) {
    char buffer[32];
    meshioplusplus::detail::snprintf_c(buffer, sizeof(buffer), "%.*f", Decimals, Value);
    return buffer;
}

std::string short_number(double Value) {
    char buffer[32];
    meshioplusplus::detail::snprintf_c(buffer, sizeof(buffer), "%.4g", Value);
    return buffer;
}

const char* edges_name(RenderEdges Edges) {
    return Edges == RenderEdges::All ? "all" : Edges == RenderEdges::Feature ? "feature" : "off";
}

// The text that keeps a status or an error short: drop the library's prefix.
std::string plain_message(std::string Text) {
    for (const char* prefix : {"meshio++: render: ", "meshio++: tui: ", "meshio++: session: ",
                               "meshio++: "}) {
        const std::string p = prefix;
        if (Text.rfind(p, 0) == 0) {
            Text.erase(0, p.size());
            break;
        }
    }
    return Text;
}

// Whitespace-separated words, with 'single' and "double" quotes keeping spaces.
std::vector<std::string> split_command(const std::string& rLine) {
    std::vector<std::string> out;
    std::string cur;
    bool in_word = false;
    char quote = 0;
    for (char c : rLine) {
        if (quote != 0) {
            if (c == quote)
                quote = 0;
            else
                cur.push_back(c);
            continue;
        }
        if (c == '\'' || c == '"') {
            quote = c;
            in_word = true;
            continue;
        }
        if (c == ' ' || c == '\t') {
            if (in_word) {
                out.push_back(cur);
                cur.clear();
                in_word = false;
            }
            continue;
        }
        cur.push_back(c);
        in_word = true;
    }
    if (quote != 0)
        throw std::invalid_argument("an unclosed quote");
    if (in_word)
        out.push_back(cur);
    return out;
}

// The name of a `--name=value` token.
std::string token_name(const std::string& rToken) {
    const std::size_t eq = rToken.find('=');
    return rToken.substr(2, eq == std::string::npos ? std::string::npos : eq - 2);
}

// What a scene is built from: every option except those the draw takes anew.
std::string scene_key(const RenderOptions& rOptions) {
    static const std::set<std::string> draw_time = {
        "azimuth", "elevation",    "roll",      "perspective", "fov",         "zoom",
        "pan-x",   "pan-y",        "axes",      "scale-bar",   "colorbar",    "cutaway",
        "cutaway-tint", "ambient", "one-sided", "background",  "point-radius", "supersample",
        "view"};
    std::string key;
    for (const std::string& token : RenderFlags::FromOptions(rOptions).Tokens())
        if (draw_time.count(token_name(token)) == 0)
            key += token + "\n";
    return key;
}

std::string cutaway_key(const RenderOptions& rOptions) {
    std::string key;
    for (const std::string& token : RenderFlags::FromOptions(rOptions).Tokens())
        if (token_name(token).rfind("cutaway", 0) == 0)
            key += token + "\n";
    return key;
}

const char* const kHelp[] = {
    " meshio++ tui                                                         ",
    " drag        orbit            wheel / + -    zoom                      ",
    " arrows      pan              r / Home       reset                     ",
    " 1 .. 7      iso +x -x +y -y +z -z      click   probe a cell          ",
    " p perspective   a axes   b scale bar   c colour bar   0  clear probes ",
    " e edges   s shading   i pin the probe   x y z  cut-away   , .  slide   ",
    " [ ]  step   { }  step by 10   space  play        :  command line      ",
    " ?  this help       q / Esc / Ctrl-C   quit                            ",
};

}  // namespace

// ---------------------------------------------------------------------------
// Construction and the prepared scenes
// ---------------------------------------------------------------------------

TuiSession::TuiSession(const Mesh& rMesh, TuiOptions Options) : mOptions(std::move(Options)) {
    mInitial = mOptions.mRender;
    apply_view_name(mInitial);
    mCamera = mInitial;
    mpSeries = mOptions.mpSeries;
    if (mOptions.mpCompare != nullptr && !is_cell_encoding(mOptions.mText.mEncoding))
        throw std::invalid_argument(
            "meshio++: tui: comparing two meshes needs a cell encoding, not a graphics protocol");
    if (mpSeries) {
        const std::size_t count = mpSeries->Count();
        if (count == 0)
            throw std::invalid_argument("meshio++: tui: the series has no steps");
        mStep = mOptions.mFollowMs > 0 ? count - 1 : 0;
        mOwnedA = mpSeries->Load(mStep);
        mpMeshA = &mOwnedA;
        mKnownFingerprint = mpSeries->Fingerprint();
    } else {
        mpMeshA = &rMesh;
    }
    if (!mOptions.mSessionPath.empty() && std::filesystem::exists(mOptions.mSessionPath))
        LoadSession(mOptions.mSessionPath);  // a bad file stops the start, by name
    PrepareScenes();
}

void TuiSession::PrepareScenes() {
    // The options a scene is built from. A series keeps one colour range across
    // its steps (the first one shown) unless told otherwise, so colours do not
    // jump from step to step.
    RenderOptions p = mCamera;
    if (mpSeries && mRangeLocked && mHaveLockedRange && !p.mVMin.has_value() &&
        !p.mVMax.has_value()) {
        p.mVMin = mLockedMin;
        p.mVMax = mLockedMax;
    }
    std::vector<Pane> panes(mOptions.mpCompare != nullptr ? 2 : 1);
    panes[0].mpMesh = mpMeshA;
    panes[0].mTitle = mOptions.mTitle;
    RenderOptions q = p;
    if (mOptions.mpCompare != nullptr) {
        panes[1].mTitle = mOptions.mCompareTitle;
        if (mOptions.mDiff) {
            const std::string& name = p.mColorBy;
            if (name.empty() || !mpMeshA->HasPointData(name) ||
                !mOptions.mpCompare->HasPointData(name) ||
                mpMeshA->PointData(name).Size() != mOptions.mpCompare->PointData(name).Size())
                throw std::invalid_argument(
                    "meshio++: tui: a difference needs color_by naming a point array that both "
                    "meshes have, over the same nodes");
            const NDArray& a = mpMeshA->PointData(name);
            const NDArray& b = mOptions.mpCompare->PointData(name);
            std::size_t k = 1;
            for (std::size_t d = 1; d < a.Ndim(); ++d)
                k *= a.Shape()[d];
            const std::size_t n = a.Size() / std::max<std::size_t>(k, 1);
            NDArray diff = NDArray::Uninit(DType::Float64, {n});
            double* out = diff.As<double>();
            for (std::size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (std::size_t c = 0; c < k; ++c) {
                    const double v = meshioplusplus::detail::read_double(b, i * k + c) -
                                     meshioplusplus::detail::read_double(a, i * k + c);
                    sum += v * v;
                }
                out[i] = std::sqrt(sum);
            }
            panes[1].mOwned =
                std::make_unique<Mesh>(meshioplusplus::detail::clone_mesh(*mOptions.mpCompare));
            panes[1].mOwned->AddPointData("diff:" + name, std::move(diff));
            panes[1].mpMesh = panes[1].mOwned.get();
            q.mColorBy = "diff:" + name;
            q.mComponent.reset();
            q.mVMin.reset();
            q.mVMax.reset();
            q.mCmap = "magma";
        } else {
            panes[1].mpMesh = mOptions.mpCompare;
        }
    }
    panes[0].mScene = prepare_render(*panes[0].mpMesh, p);
    if (panes.size() > 1) {
        panes[1].mScene = prepare_render(*panes[1].mpMesh, q);
        // One colour range over both, read off a small frame of each.
        if (mOptions.mSharedRange && !mOptions.mDiff && !p.mVMin.has_value() &&
            !p.mVMax.has_value()) {
            RenderOptions probe = p;
            probe.mWidth = 8;
            probe.mHeight = 8;
            const Frame fa = render_scene(panes[0].mScene, probe);
            const Frame fb = render_scene(panes[1].mScene, probe);
            if (fa.mColored && fb.mColored) {
                p.mVMin = std::min(fa.mVMin, fb.mVMin);
                p.mVMax = std::max(fa.mVMax, fb.mVMax);
                panes[0].mScene = prepare_render(*panes[0].mpMesh, p);
                panes[1].mScene = prepare_render(*panes[1].mpMesh, p);
            }
        }
    }
    mPanes = std::move(panes);
    mPrepared = mCamera;
    mPreparedSmooth = mCamera.mShading == RenderShading::Smooth;
    mReprepare = false;
}

void TuiSession::RollBackPreparation() {
    // A change the scene could not be built with (a colormap that does not exist,
    // a field the mesh lacks) must not stay in the options: put back what the
    // scene was prepared with and keep the camera and the draw-time toggles the
    // user moved meanwhile.
    const RenderOptions bad = mCamera;
    mCamera = mPrepared;
    mCamera.mAzimuth = bad.mAzimuth;
    mCamera.mElevation = bad.mElevation;
    mCamera.mRoll = bad.mRoll;
    mCamera.mProjection = bad.mProjection;
    mCamera.mFovDeg = bad.mFovDeg;
    mCamera.mZoom = bad.mZoom;
    mCamera.mPanX = bad.mPanX;
    mCamera.mPanY = bad.mPanY;
    mCamera.mAxes = bad.mAxes;
    mCamera.mScaleBar = bad.mScaleBar;
    mCamera.mColorbar = bad.mColorbar;
    mCamera.mCutaways = bad.mCutaways;
    mCamera.mCutawayTint = bad.mCutawayTint;
    mCamera.mAmbient = bad.mAmbient;
    mCamera.mTwoSided = bad.mTwoSided;
    mCamera.mBackground = bad.mBackground;
    mCamera.mPointRadius = bad.mPointRadius;
    mReprepare = false;
}

bool TuiSession::LoadStep(std::size_t Index) {
    if (!mpSeries || Index >= mpSeries->Count())
        return false;
    try {
        mOwnedA = mpSeries->Load(Index);
    } catch (const std::exception& rErr) {
        mMessage = "step " + std::to_string(Index + 1) + ": " + plain_message(rErr.what());
        return false;
    }
    mpMeshA = &mOwnedA;
    mStep = Index;
    mReprepare = true;
    return true;
}

bool TuiSession::StepBy(long long Delta) {
    if (!mpSeries) {
        mMessage = "one mesh, no steps (give a glob or several files)";
        return true;
    }
    const long long count = static_cast<long long>(mpSeries->Count());
    const long long target =
        std::clamp<long long>(static_cast<long long>(mStep) + Delta, 0, count - 1);
    if (target == static_cast<long long>(mStep))
        return false;
    return LoadStep(static_cast<std::size_t>(target));
}

// ---------------------------------------------------------------------------
// The camera and the keys
// ---------------------------------------------------------------------------

bool TuiSession::ResetCamera() {
    // The camera returns to where it started; the field, edges and shading keep
    // what the user set, since the prepared scene was built with them. The
    // cut-aways go.
    const RenderOptions keep = mCamera;
    mCamera = mInitial;
    mCamera.mEdges = keep.mEdges;
    mCamera.mShading = keep.mShading;
    mCamera.mCutaways.clear();
    mCutSign = {0, 0, 0};
    mCutActive = -1;
    return true;
}

void TuiSession::SyncCutaways() {
    mCamera.mCutaways.clear();
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (mCutSign[axis] == 0)
            continue;
        RenderCutaway plane;
        plane.mPoint = {0.0, 0.0, 0.0};
        plane.mNormal = {0.0, 0.0, 0.0};
        plane.mPoint[axis] = mCutAt[axis];
        plane.mNormal[axis] = static_cast<double>(mCutSign[axis]);
        mCamera.mCutaways.push_back(plane);
    }
}

bool TuiSession::ToggleCutaway(int Axis) {
    const std::size_t a = static_cast<std::size_t>(Axis);
    // off -> keeps the + side -> keeps the - side -> off
    const int next = mCutSign[a] == 0 ? 1 : mCutSign[a] == 1 ? -1 : 0;
    int active = 0;
    for (int s : mCutSign)
        active += s != 0 ? 1 : 0;
    if (mCutSign[a] == 0 && active >= 2) {
        mMessage = "at most two cut-away planes: clear one first";
        return true;
    }
    if (mCutSign[a] == 0 && !mPanes.empty()) {
        std::array<double, 3> centre{};
        double radius = 0.0;
        if (render_scene_bounds(mPanes[0].mScene, centre, radius))
            mCutAt[a] = centre[a];
    }
    mCutSign[a] = next;
    mCutActive = next == 0 ? -1 : Axis;
    SyncCutaways();
    static const char axes[] = {'x', 'y', 'z'};
    mMessage = next == 0 ? std::string("cut-away ") + axes[a] + " off"
                         : std::string("cut-away: keeps ") + axes[a] + (next > 0 ? " >= " : " <= ") +
                               short_number(mCutAt[a]);
    return true;
}

bool TuiSession::MoveCutaway(double Direction) {
    if (mCutActive < 0 || mPanes.empty())
        return false;
    std::array<double, 3> centre{};
    double radius = 0.0;
    if (!render_scene_bounds(mPanes[0].mScene, centre, radius))
        return false;
    const std::size_t a = static_cast<std::size_t>(mCutActive);
    mCutAt[a] += Direction * radius / 25.0;
    SyncCutaways();
    static const char axes[] = {'x', 'y', 'z'};
    mMessage = std::string("cut-away ") + axes[a] + (mCutSign[a] > 0 ? " >= " : " <= ") +
               short_number(mCutAt[a]);
    return true;
}

bool TuiSession::Handle(const InputEvent& rEvent, const TerminalSize& rSize) {
    const int cols = std::max(1, rSize.mCols);
    const int pic_rows = std::max(1, rSize.mRows - 1 - mNoteRows - mProbeRows);
    auto zoom_by = [&](double Factor) {
        mCamera.mZoom = std::clamp(mCamera.mZoom * Factor, kMinZoom, kMaxZoom);
        return true;
    };
    auto orbit = [&](int Dx, int Dy) {
        mCamera.mAzimuth -= Dx * 360.0 / cols;
        mCamera.mAzimuth = std::remainder(mCamera.mAzimuth, 360.0);
        mCamera.mElevation = std::clamp(mCamera.mElevation + Dy * 180.0 / pic_rows, -90.0, 90.0);
        return Dx != 0 || Dy != 0;
    };
    auto pan = [&](double Dx, double Dy) {
        mCamera.mPanX += Dx;
        mCamera.mPanY += Dy;
        return Dx != 0.0 || Dy != 0.0;
    };

    if (rEvent.mKind == EventKind::Mouse) {
        switch (rEvent.mAction) {
            case MouseAction::Press:
                mDragging = true;
                mDragMoved = false;
                mDragButton = rEvent.mButton;
                mLastX = rEvent.mX;
                mLastY = rEvent.mY;
                return false;
            case MouseAction::Release: {
                const bool click = mDragging && !mDragMoved && mDragButton == 0;
                mDragging = false;
                if (click) {
                    ProbeAt(rEvent.mX, rEvent.mY);
                    return true;
                }
                return false;
            }
            case MouseAction::Drag: {
                if (!mDragging) {
                    mDragging = true;
                    mDragMoved = true;
                    mDragButton = rEvent.mButton;
                    mLastX = rEvent.mX;
                    mLastY = rEvent.mY;
                    return false;
                }
                const int dx = rEvent.mX - mLastX;
                const int dy = rEvent.mY - mLastY;
                mLastX = rEvent.mX;
                mLastY = rEvent.mY;
                if (dx != 0 || dy != 0)
                    mDragMoved = true;
                if (mDragButton == 0 && !rEvent.mShift)
                    return orbit(dx, dy);
                return pan(static_cast<double>(dx) / cols, -static_cast<double>(dy) / pic_rows);
            }
            case MouseAction::WheelUp:
                return zoom_by(kZoomStep);
            case MouseAction::WheelDown:
                return zoom_by(1.0 / kZoomStep);
            case MouseAction::Move:
                return false;
        }
        return false;
    }
    if (rEvent.mKind == EventKind::Paste) {
        if (mPrompt) {
            for (char c : rEvent.mText)
                if (c != '\n' && c != '\r')
                    mPromptText.push_back(c);
            return true;
        }
        return false;
    }

    if (mPrompt)
        return HandlePromptKey(rEvent);

    mMessage.clear();
    if (rEvent.mKey == KeyCode::Escape || rEvent.IsChar('q') || rEvent.IsCtrl('c') ||
        rEvent.IsCtrl('d')) {
        mQuit = true;
        return false;
    }
    switch (rEvent.mKey) {
        case KeyCode::Up:
            return pan(0.0, kPanStep);
        case KeyCode::Down:
            return pan(0.0, -kPanStep);
        case KeyCode::Left:
            return pan(-kPanStep, 0.0);
        case KeyCode::Right:
            return pan(kPanStep, 0.0);
        case KeyCode::PageUp:
            return zoom_by(kZoomStep);
        case KeyCode::PageDown:
            return zoom_by(1.0 / kZoomStep);
        case KeyCode::Home:
            return ResetCamera();
        default:
            break;
    }
    if (rEvent.mKey != KeyCode::Char || rEvent.mCtrl || rEvent.mAlt)
        return false;
    const std::uint32_t ch = rEvent.mChar;
    if (ch >= '1' && ch <= '7') {
        const Preset& preset = kPresets[ch - '1'];
        mCamera.mAzimuth = preset.mAzimuth;
        mCamera.mElevation = preset.mElevation;
        mCamera.mRoll = 0.0;
        mMessage = std::string("view ") + preset.mName;
        return true;
    }
    switch (ch) {
        case '+':
        case '=':
            return zoom_by(kZoomStep);
        case '-':
        case '_':
            return zoom_by(1.0 / kZoomStep);
        case 'r':
            return ResetCamera();
        case 'p':
            mCamera.mProjection = mCamera.mProjection == RenderProjection::Perspective
                                      ? RenderProjection::Orthographic
                                      : RenderProjection::Perspective;
            return true;
        case 'a':
            mCamera.mAxes = !mCamera.mAxes;
            return true;
        case 'b':
            mCamera.mScaleBar = !mCamera.mScaleBar;
            return true;
        case 'c':
            mCamera.mColorbar = !mCamera.mColorbar;
            return true;
        case 'e':
            mCamera.mEdges = mCamera.mEdges == RenderEdges::None  ? RenderEdges::All
                             : mCamera.mEdges == RenderEdges::All ? RenderEdges::Feature
                                                                  : RenderEdges::None;
            mReprepare = true;  // edge lines are part of the prepared scene
            mMessage = std::string("edges ") + edges_name(mCamera.mEdges);
            return true;
        case 's':
            mCamera.mShading = mCamera.mShading == RenderShading::Flat     ? RenderShading::Smooth
                               : mCamera.mShading == RenderShading::Smooth ? RenderShading::None
                                                                           : RenderShading::Flat;
            if (mCamera.mShading == RenderShading::Smooth && !mPreparedSmooth)
                mReprepare = true;  // smooth normals are computed when preparing
            mMessage = "shading " + std::string(mCamera.mShading == RenderShading::Smooth ? "smooth"
                                                : mCamera.mShading == RenderShading::Flat ? "flat"
                                                                                          : "none");
            return true;
        case 'x':
            return ToggleCutaway(0);
        case 'y':
            return ToggleCutaway(1);
        case 'z':
            return ToggleCutaway(2);
        case ',':
            return MoveCutaway(-1.0);
        case '.':
            return MoveCutaway(1.0);
        case 'i':
            if (mProbe.mCell < 0) {
                mMessage = "click a cell first";
                return true;
            }
            if (mPins.size() == 2)
                mPins.erase(mPins.begin());
            mPins.push_back(mProbe);
            RebuildProbeLines();
            mMessage = "pinned";
            return true;
        case '0':
            mProbe = Probe();
            mPins.clear();
            mCutSign = {0, 0, 0};
            mCutActive = -1;
            SyncCutaways();
            RebuildProbeLines();
            mMessage = "cleared probes and cut-aways";
            return true;
        case '[':
            return StepBy(-1);
        case ']':
            return StepBy(1);
        case '{':
            return StepBy(-10);
        case '}':
            return StepBy(10);
        case ' ':
            if (!mpSeries || mpSeries->Count() < 2) {
                mMessage = "one step: nothing to play";
                return true;
            }
            mPlaying = !mPlaying;
            mNextPlayMs = mNowMs;
            mMessage = mPlaying ? "playing" : "paused";
            return true;
        case ':':
            mPrompt = true;
            mPromptText.clear();
            return true;
        case '?':
        case 'h':
            mHelp = !mHelp;
            return true;
        default:
            return false;
    }
}

bool TuiSession::HandlePromptKey(const InputEvent& rEvent) {
    if (rEvent.mKey == KeyCode::Escape || rEvent.IsCtrl('c')) {
        mPrompt = false;
        mPromptText.clear();
        return true;
    }
    if (rEvent.mKey == KeyCode::Enter) {
        const std::string line = mPromptText;
        mPrompt = false;
        mPromptText.clear();
        mMessage = Execute(line);
        return true;
    }
    if (rEvent.mKey == KeyCode::Backspace) {
        if (!mPromptText.empty()) {
            // Remove one UTF-8 code point.
            std::size_t i = mPromptText.size() - 1;
            while (i > 0 && (static_cast<unsigned char>(mPromptText[i]) & 0xC0) == 0x80)
                --i;
            mPromptText.erase(i);
        }
        return true;
    }
    if (rEvent.IsCtrl('u')) {
        mPromptText.clear();
        return true;
    }
    if (rEvent.mKey == KeyCode::Char && !rEvent.mCtrl && !rEvent.mAlt && rEvent.mChar >= 0x20) {
        // UTF-8 encode the code point.
        const std::uint32_t cp = rEvent.mChar;
        if (cp < 0x80) {
            mPromptText.push_back(static_cast<char>(cp));
        } else if (cp < 0x800) {
            mPromptText.push_back(static_cast<char>(0xC0 | (cp >> 6)));
            mPromptText.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else if (cp < 0x10000) {
            mPromptText.push_back(static_cast<char>(0xE0 | (cp >> 12)));
            mPromptText.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            mPromptText.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else {
            mPromptText.push_back(static_cast<char>(0xF0 | (cp >> 18)));
            mPromptText.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
            mPromptText.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            mPromptText.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Probing
// ---------------------------------------------------------------------------

void TuiSession::ProbeAt(int X, int Y) {
    for (const Pane& pane : mPanes) {
        if (X - 1 < pane.mColOffset || X - 1 >= pane.mColOffset + pane.mCols)
            continue;
        const int col = X - 1 - pane.mColOffset;
        const int row = Y - 1;
        const int width = pane.mFrame.mWidth;
        const int height = pane.mFrame.mHeight;
        if (width <= 0 || row < 0 || row * pane.mGy >= height) {
            mMessage = "nothing there";
            return;
        }
        std::int64_t id = -1;
        for (int sy = 0; sy < pane.mGy && id < 0; ++sy)
            for (int sx = 0; sx < pane.mGx && id < 0; ++sx) {
                const int px = col * pane.mGx + sx;
                const int py = row * pane.mGy + sy;
                if (px < width && py < height)
                    id = pane.mFrame.mCellIds[static_cast<std::size_t>(py) *
                                                  static_cast<std::size_t>(width) +
                                              static_cast<std::size_t>(px)];
            }
        if (id < 0) {
            mMessage = "nothing there";
            return;
        }
        mProbe = probe_cell(*pane.mpMesh, id);
        mMessage.clear();
        RebuildProbeLines();
        return;
    }
}

void TuiSession::RebuildProbeLines() {
    // This probe (up to three lines), then the pins and their difference.
    mProbeLines.clear();
    for (std::size_t i = 0; i < mProbe.mLines.size() && i < 3; ++i)
        mProbeLines.push_back(mProbe.mLines[i]);
    for (std::size_t k = 0; k < mPins.size(); ++k) {
        std::string line = std::string("pin ") + static_cast<char>('A' + k) + ": cell " +
                           std::to_string(mPins[k].mCell);
        int shown = 0;
        for (const auto& [name, value] : mPins[k].mValues) {
            if (shown++ >= 4)
                break;
            line += "  " + name + "=" + short_number(value);
        }
        mProbeLines.push_back(line);
    }
    if (mPins.size() == 2) {
        const std::string diff = probe_difference(mPins[0], mPins[1]);
        mProbeLines.push_back("B - A: " +
                              (diff.empty() ? std::string("(no field in common)") : diff));
    }
    const int rows = std::min<int>(static_cast<int>(mProbeLines.size()), kMaxProbeRows);
    if (rows != mProbeRows) {
        mProbeRows = rows;
        mFull = true;  // the picture changes height
    }
}

// ---------------------------------------------------------------------------
// Commands, sessions, snapshots
// ---------------------------------------------------------------------------

void TuiSession::ApplyFlags(RenderFlags& rFlags, bool Prepare) {
    RenderOptions next = rFlags.ToOptions();
    apply_view_name(next);
    const bool scene_changed = scene_key(next) != scene_key(mCamera);
    const bool cut_changed = cutaway_key(next) != cutaway_key(mCamera);
    const RenderOptions before = mCamera;
    const std::array<int, 3> cut_sign = mCutSign;
    const int cut_active = mCutActive;
    if (cut_changed) {
        mCutSign = {0, 0, 0};
        mCutActive = -1;
    }
    mCamera = next;
    if (!scene_changed)
        return;
    if (!Prepare) {
        mReprepare = true;
        return;
    }
    // A change the scene is built from is built now, so one that cannot be (a
    // colormap that does not exist, an array the mesh lacks) is refused with its
    // reason and leaves everything before it as it was.
    if (mpIo != nullptr && !mPanes.empty()) {
        TerminalSize size;
        if (mpIo->Size(size))
            mpIo->Write(at_row(size.mRows) + "\x1b[7m" + fit_columns(" preparing...", size.mCols) +
                        "\x1b[0m");
    }
    try {
        PrepareScenes();
    } catch (...) {
        mCamera = before;
        mCutSign = cut_sign;
        mCutActive = cut_active;
        throw;
    }
}

std::string TuiSession::WriteSnapshot(const std::string& rPath) {
    RenderOptions o = mCamera;
    o.mWidth = 800;
    o.mHeight = 600;
    TextOptions t = mOptions.mText;
    t.mCols = 100;
    t.mRows = 40;
    t.mNotes = true;
    meshioplusplus::write_snapshot(rPath, *mpMeshA, o, t, SnapshotOptions{});
    return "wrote " + rPath;
}

void TuiSession::SaveSession(const std::string& rPath) {
    SessionFile session;
    session.mFlags = RenderFlags::FromOptions(mCamera).Tokens();
    session.mStep = static_cast<long long>(mStep);
    const std::string text = session_to_json(session);
    auto out = meshioplusplus::detail::make_classic_ofstream(
        rPath, std::ios_base::out | std::ios_base::binary);
    if (!out)
        throw std::runtime_error("cannot write the session file '" + rPath + "'");
    out << text;
    out.flush();
    if (!out)
        throw std::runtime_error("cannot write the session file '" + rPath + "'");
}

void TuiSession::LoadSession(const std::string& rPath) {
    auto in = meshioplusplus::detail::make_classic_ifstream(
        rPath, std::ios_base::in | std::ios_base::binary);
    if (!in)
        throw std::runtime_error("cannot read the session file '" + rPath + "'");
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const SessionFile session = session_from_json(text);
    RenderFlags flags = RenderFlags::FromTokens(session.mFlags);
    ApplyFlags(flags, /*Prepare=*/false);  // the constructor prepares once, after this
    if (mpSeries && session.mStep > 0)
        LoadStep(std::min<std::size_t>(static_cast<std::size_t>(session.mStep),
                                       mpSeries->Count() - 1));
}

std::string TuiSession::Execute(const std::string& rLine) {
    try {
        const std::vector<std::string> t = split_command(rLine);
        if (t.empty())
            return "";
        std::string name = t[0];
        if (!name.empty() && name[0] == ':')
            name.erase(0, 1);
        auto need = [&](std::size_t n, const std::string& rUsage) {
            if (t.size() < n + 1)
                throw std::invalid_argument("usage: :" + rUsage);
        };

        if (name == "q" || name == "quit") {
            mQuit = true;
            return "";
        }
        if (name == "help" || name == "?") {
            mHelp = !mHelp;
            return "";
        }
        if (name == "reset") {
            ResetCamera();
            return "reset";
        }
        if (name == "w" || name == "write") {
            need(1, "w FILE   (.png .txt .ansi .html .cast)");
            return WriteSnapshot(t[1]);
        }
        if (name == "session") {
            need(2, "session save|load FILE");
            if (t[1] == "save") {
                SaveSession(t[2]);
                return "saved " + t[2];
            }
            if (t[1] == "load") {
                LoadSession(t[2]);
                return "loaded " + t[2];
            }
            throw std::invalid_argument("usage: :session save|load FILE");
        }
        if (name == "step") {
            need(1, "step N | next | prev | first | last");
            if (!mpSeries)
                return "one mesh, no steps (give a glob or several files)";
            const long long count = static_cast<long long>(mpSeries->Count());
            long long target = static_cast<long long>(mStep);
            if (t[1] == "next")
                target = std::min(count - 1, target + 1);
            else if (t[1] == "prev")
                target = std::max(0LL, target - 1);
            else if (t[1] == "first")
                target = 0;
            else if (t[1] == "last")
                target = count - 1;
            else {
                const auto v = cli_parse_numbers(t[1], "step");
                if (v.size() != 1 || !v[0].has_value())
                    throw std::invalid_argument(
                        "usage: :step N (1-based) | next | prev | first | last");
                target = std::clamp(static_cast<long long>(*v[0]) - 1, 0LL, count - 1);
            }
            if (target != static_cast<long long>(mStep))
                LoadStep(static_cast<std::size_t>(target));
            return mMessage;
        }
        if (name == "probe") {
            if (t.size() > 1 && t[1] == "off") {
                mProbe = Probe();
                mPins.clear();
                RebuildProbeLines();
                return "probe off";
            }
            throw std::invalid_argument("usage: :probe off   (click a cell to probe it)");
        }
        if (name == "pin") {
            if (mProbe.mCell < 0)
                return "click a cell first";
            if (mPins.size() == 2)
                mPins.erase(mPins.begin());
            mPins.push_back(mProbe);
            RebuildProbeLines();
            return "pinned";
        }

        // Everything below changes the render flags: take the state as flags,
        // edit it, and derive the options again, so a command is exactly a
        // `snapshot` flag.
        RenderFlags flags = RenderFlags::FromOptions(mCamera);
        if (name == "unset") {
            need(1, "unset FLAG");
            if (!RenderFlags::Known(t[1]))
                throw std::invalid_argument("unknown flag '" + t[1] + "'");
            flags.Unset(t[1]);
        } else if (name == "range") {
            need(1, "range LOW HIGH | auto | lock");
            if (t[1] == "auto" || t[1] == "lock") {
                flags.Unset("vmin");
                flags.Unset("vmax");
                mRangeLocked = t[1] == "lock";
                mHaveLockedRange = false;
                mReprepare = true;
            } else {
                need(2, "range LOW HIGH | auto | lock");
                flags.Set("vmin", t[1]);
                flags.Set("vmax", t[2]);
            }
        } else if (name == "clip") {
            need(1, "clip AXIS:OFFSET | off     (AXIS is +x -x +y -y +z -z; `:clip +x 0.5` too)");
            if (t[1] == "off") {
                flags.Unset("cutaway");
                flags.Unset("cutaway-tint");
            } else {
                const std::string text = t.size() > 2 ? t[1] + ":" + t[2] : t[1];
                std::vector<RenderCutaway> planes = mCamera.mCutaways;
                planes.push_back(cli_parse_cutaway(text));
                if (planes.size() > 2)
                    planes.erase(planes.begin());  // a third replaces the first
                RenderOptions carrier;
                carrier.mCutaways = planes;
                const RenderFlags planes_as_flags = RenderFlags::FromOptions(carrier);
                flags.Unset("cutaway");
                for (const std::string& value : planes_as_flags.Values("cutaway"))
                    flags.Add("cutaway", value);
            }
        } else if (name == "color" || name == "color-by") {
            need(1, "color NAME | none");
            if (t[1] == "none")
                flags.Unset("color-by");
            else
                flags.Set("color-by", t[1]);
        } else if (name == "view") {
            need(1, "view iso|+x|-x|+y|-y|+z|-z");
            flags.Unset("azimuth");
            flags.Unset("elevation");
            flags.Set("view", t[1]);
        } else if (RenderFlags::Known(name)) {
            if (RenderFlags::TakesValue(name)) {
                need(1, name + " VALUE");
                const bool clearable = name == "color-by" || name == "expr" || name == "vectors" ||
                                       name == "streamlines" || name == "warp" || name == "reduce" ||
                                       name == "quality-metric";
                if (t[1] == "none" && clearable)
                    flags.Unset(name);
                else
                    flags.Set(name, t[1]);
            } else {
                const bool on = t.size() > 1 ? (t[1] == "on" || t[1] == "true" || t[1] == "1")
                                             : !flags.Has(name);
                flags.SetFlag(name, on);
            }
        } else {
            return "unknown command '" + name +
                   "' (:help lists the keys; every snapshot flag works)";
        }
        ApplyFlags(flags);
        return "";
    } catch (const std::exception& rErr) {
        return plain_message(rErr.what());
    }
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

void TuiSession::Draw(TuiIo& rIo, const TerminalSize& rSize, bool Full) {
    const int cols = rSize.mCols;
    const int rows = rSize.mRows;
    std::string out;
    if (Full)
        out += "\x1b[0m\x1b[2J";
    if (mReprepare) {
        // The preparation can take seconds on a large volume: say so first.
        rIo.Write(at_row(rows) + "\x1b[7m" + fit_columns(" preparing...", cols) + "\x1b[0m");
        try {
            PrepareScenes();
        } catch (const std::exception& rErr) {
            RollBackPreparation();
            mMessage = plain_message(rErr.what());
        }
    }
    if (cols < 8 || rows < 4) {
        out += "\x1b[2J" + at_row(1) + "terminal too small";
        rIo.Write(out);
        mPrevious = TextGrid{};
        return;
    }

    const bool cells = is_cell_encoding(mOptions.mText.mEncoding);
    const bool compare = mPanes.size() > 1;
    std::string graphics;
    int pic_rows = 1;
    const int probe_rows = std::min(mProbeRows, std::max(0, rows - 6));
    for (int pass = 0; pass < 2; ++pass) {
        pic_rows = std::max(1, rows - 1 - mNoteRows - probe_rows);
        TextOptions text = mOptions.mText;
        text.mRows = pic_rows;
        if (mOptions.mCellAspectFromTerminal && rSize.mPixelWidth > 0 && rSize.mPixelHeight > 0) {
            text.mCellPixelWidth = std::max(1, rSize.mPixelWidth / cols);
            text.mCellPixelHeight = std::max(1, rSize.mPixelHeight / rows);
            text.mCellAspect = static_cast<double>(rSize.mPixelHeight * cols) /
                               static_cast<double>(rSize.mPixelWidth * rows);
        }
        const std::array<int, 2> grid_px = cell_pixels(text.mEncoding);
        const int left = compare ? (cols - 1) / 2 : cols;
        bool failed = false;
        for (std::size_t k = 0; k < mPanes.size(); ++k) {
            Pane& pane = mPanes[k];
            pane.mCols = k == 0 ? left : cols - 1 - left;
            pane.mColOffset = k == 0 ? 0 : left + 1;
            pane.mGx = grid_px[0];
            pane.mGy = grid_px[1];
            TextOptions pane_text = text;
            pane_text.mCols = pane.mCols;
            double aspect = 1.0;
            const std::array<int, 2> size = text_frame_size(pane_text, aspect);
            RenderOptions r = mCamera;
            r.mWidth = size[0];
            r.mHeight = size[1];
            r.mPixelAspect = aspect;
            try {
                pane.mFrame = render_scene(pane.mScene, r, mOptions.mFixedFit);
                if (cells) {
                    pane.mGrid = encode_cells(pane.mFrame, pane_text);
                } else {
                    TextOptions plain = pane_text;
                    plain.mNotes = false;
                    graphics = encode_text(pane.mFrame, plain);
                }
            } catch (const std::exception& rErr) {
                mMessage = plain_message(rErr.what());
                failed = true;
                break;
            }
        }
        if (failed)
            break;
        // The frame's notes go under the picture; if their number changed, fit
        // the picture to what is left and draw once more.
        const int wanted = static_cast<int>(mPanes[0].mFrame.mNotes.size());
        if (wanted == mNoteRows || pass == 1)
            break;
        mNoteRows = std::min(wanted, std::max(0, rows - 4 - probe_rows));
        Full = true;
        out = "\x1b[0m\x1b[2J";
    }

    // A series keeps the colour range of the first step it shows.
    if (mpSeries && mRangeLocked && !mHaveLockedRange && mPanes[0].mFrame.mColored &&
        !mCamera.mVMin.has_value() && !mCamera.mVMax.has_value()) {
        mHaveLockedRange = true;
        mLockedMin = mPanes[0].mFrame.mVMin;
        mLockedMax = mPanes[0].mFrame.mVMax;
    }

    TextGrid grid;
    if (cells) {
        if (!compare) {
            grid = mPanes[0].mGrid;
        } else {
            grid.mCols = cols;
            grid.mRows = mPanes[0].mGrid.mRows;
            grid.mCells.assign(static_cast<std::size_t>(grid.mCols) *
                                   static_cast<std::size_t>(grid.mRows),
                               TextCell{});
            for (int r = 0; r < grid.mRows; ++r) {
                for (std::size_t k = 0; k < mPanes.size(); ++k) {
                    const TextGrid& g = mPanes[k].mGrid;
                    if (g.mRows != grid.mRows)
                        continue;
                    for (int c = 0; c < g.mCols; ++c)
                        grid.mCells[static_cast<std::size_t>(r) * static_cast<std::size_t>(cols) +
                                    static_cast<std::size_t>(mPanes[k].mColOffset + c)] =
                            g.mCells[static_cast<std::size_t>(r) *
                                         static_cast<std::size_t>(g.mCols) +
                                     static_cast<std::size_t>(c)];
                }
                TextCell bar;
                bar.mGlyph = 0x2502;  // the line between the two meshes
                bar.mFgSet = true;
                bar.mFg = {110, 110, 110};
                grid.mCells[static_cast<std::size_t>(r) * static_cast<std::size_t>(cols) +
                            static_cast<std::size_t>(mPanes[0].mCols)] = bar;
            }
        }
        // After a clear the screen is blank, so a full repaint only writes what is
        // not.
        TextGrid blank;
        if (Full) {
            blank.mCols = grid.mCols;
            blank.mRows = grid.mRows;
            blank.mCells.assign(grid.mCells.size(), TextCell{});
        }
        out += encode_cells_update(Full ? blank : mPrevious, grid, mOptions.mText.mDepth, 1, 1);
        mPrevious = grid;
        mReport.mScreen = grid;
    } else {
        out += "\x1b[H";
        if (mOptions.mText.mEncoding == TextEncoding::Kitty)
            out += "\x1b_Ga=d\x1b\\";  // the previous image goes before the new one
        out += graphics;
    }
    if (mHelp) {
        const int lines = static_cast<int>(sizeof(kHelp) / sizeof(kHelp[0]));
        for (int i = 0; i < lines && i < pic_rows; ++i)
            out += at_row(1 + i) + "\x1b[7m" + fit_columns(kHelp[i], std::min(cols, 72)) + "\x1b[0m";
        mHelpWasShown = true;
    } else if (mHelpWasShown) {
        mHelpWasShown = false;
        mPrevious = TextGrid{};  // the help covered cells the diff thinks are intact
        if (cells) {
            TextGrid blank;
            blank.mCols = grid.mCols;
            blank.mRows = grid.mRows;
            blank.mCells.assign(grid.mCells.size(), TextCell{});
            out = "\x1b[0m\x1b[2J" + encode_cells_update(blank, grid, mOptions.mText.mDepth);
            mPrevious = grid;
        }
    }

    // The titles of a comparison, the notes, the probe lines, then the status.
    std::vector<std::string> lines;
    std::vector<std::string> notes = mPanes[0].mFrame.mNotes;
    if (compare && !notes.empty()) {
        const std::string titles = (mPanes[0].mTitle.empty() ? "A" : mPanes[0].mTitle) + "  |  " +
                                   (mPanes[1].mTitle.empty() ? "B" : mPanes[1].mTitle);
        notes[0] = titles + "   " + notes[0];
    }
    for (int i = 0; i < mNoteRows && i < static_cast<int>(notes.size()); ++i)
        lines.push_back(fit_columns(notes[static_cast<std::size_t>(i)], cols));
    while (static_cast<int>(lines.size()) < mNoteRows)
        lines.push_back(fit_columns("", cols));
    for (int i = 0; i < probe_rows; ++i)
        lines.push_back(fit_columns(i < static_cast<int>(mProbeLines.size())
                                        ? mProbeLines[static_cast<std::size_t>(i)]
                                        : std::string(),
                                    cols));
    std::string status;
    if (mPrompt) {
        status = " :" + mPromptText + "_";
    } else {
        status = " " + (mOptions.mTitle.empty() ? std::string("meshio++") : mOptions.mTitle);
        if (mpSeries)
            status += "  " + mpSeries->Label(mStep) + (mPlaying ? " >" : "");
        status += "  az " + number(mCamera.mAzimuth, 0) + "  el " + number(mCamera.mElevation, 0) +
                  "  zoom " + number(mCamera.mZoom, 2) + "  " +
                  (mCamera.mProjection == RenderProjection::Perspective ? "persp" : "ortho");
        if (!mCamera.mCutaways.empty())
            status += "  cut " + std::to_string(mCamera.mCutaways.size());
        if (!mMessage.empty())
            status += "  | " + mMessage;
        status += "  | ? help  q quit";
    }
    mReport.mStatus = status;
    lines.push_back("\x1b[7m" + fit_columns(status, cols) + "\x1b[0m");
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const bool is_status = i + 1 == lines.size();
        if (!Full && i < mPreviousLines.size() && mPreviousLines[i] == lines[i] && !is_status)
            continue;
        const int row = is_status ? rows : pic_rows + 1 + static_cast<int>(i);
        out += at_row(row) + lines[i] + "\x1b[0m";
    }
    mPreviousLines = lines;
    rIo.Write(out);
    ++mReport.mFrames;
}

// ---------------------------------------------------------------------------
// The loop
// ---------------------------------------------------------------------------

void TuiSession::Tick() {
    // Playing: the next step every 1/fps seconds, stopping at the last.
    if (mPlaying && mpSeries && mNowMs >= mNextPlayMs) {
        const double fps = mOptions.mPlayFps > 0.0 ? mOptions.mPlayFps : 4.0;
        mNextPlayMs = mNowMs + static_cast<long long>(1000.0 / fps);
        if (!StepBy(1))
            mPlaying = false;
    }
    // Following: notice a new or changing last step, wait for it to settle, read
    // it, and keep what was on screen when the read fails (a half-written file is
    // expected, not fatal).
    if (mOptions.mFollowMs > 0 && mpSeries) {
        if (!mSettling && mNowMs >= mNextFollowMs) {
            mNextFollowMs = mNowMs + mOptions.mFollowMs;
            const std::string fp = mpSeries->Fingerprint();
            if (fp != mKnownFingerprint) {
                mSettling = true;
                mSettlingFingerprint = fp;
                mSettleUntilMs = mNowMs + mOptions.mFollowSettleMs;
            }
        } else if (mSettling && mNowMs >= mSettleUntilMs) {
            const std::string fp = mpSeries->Fingerprint();
            if (fp != mSettlingFingerprint) {
                mSettlingFingerprint = fp;  // still being written
                mSettleUntilMs = mNowMs + mOptions.mFollowSettleMs;
            } else {
                mSettling = false;
                mpSeries->Refresh();
                const std::size_t count = mpSeries->Count();
                if (count > 0 && LoadStep(count - 1)) {
                    mKnownFingerprint = fp;
                    mMessage = "followed: " + mpSeries->Label(mStep);
                } else {
                    mMessage = "waiting for a complete step";
                }
            }
        }
    }
}

TuiReport TuiSession::Run(TuiIo& rIo) {
    mpIo = &rIo;
    TerminalSize size;
    if (!rIo.Size(size)) {
        size.mCols = 80;
        size.mRows = 24;
    }
    InputParser parser;
    bool dirty = true;
    bool full = true;
    while (true) {
        if (const int sig = rIo.TerminationSignal(); sig != 0) {
            mReport.mExit = 128 + sig;
            break;
        }
        if (rIo.TakeResize()) {
            rIo.Size(size);
            dirty = full = true;
        }
        mNowMs = rIo.NowMs();
        const std::size_t step_before = mStep;
        const std::string message_before = mMessage;
        const bool playing_before = mPlaying;
        Tick();
        if (mStep != step_before || mMessage != message_before || mPlaying != playing_before)
            dirty = true;
        if (mFull) {
            full = true;
            dirty = true;
            mFull = false;
        }
        if (dirty || mReprepare) {
            Draw(rIo, size, full);
            dirty = full = false;
        }
        std::string bytes;
        const bool alive = rIo.Read(bytes, mOptions.mPollMs);
        parser.Feed(bytes);
        for (const InputEvent& event : parser.Take(bytes.empty())) {
            if (Handle(event, size))
                dirty = true;
            if (mQuit)
                break;
        }
        if (mQuit || !alive)
            break;
    }
    // Keys that arrived together with the one that ended the session (a recorded
    // script, a pasted burst) still count: show where they left the camera.
    if ((dirty || mFull || mReprepare) && mReport.mExit == 0)
        Draw(rIo, size, full || mFull);
    mpIo = nullptr;
    mReport.mFinal = mCamera;
    mReport.mStep = mStep;
    mReport.mProbe = mProbeLines;
    if (!mOptions.mSessionPath.empty()) {
        try {
            SaveSession(mOptions.mSessionPath);
        } catch (const std::exception& rErr) {
            mReport.mError = plain_message(rErr.what());
        }
    }
    return mReport;
}

TuiReport run_on_terminal(const Mesh& rMesh, TuiOptions Options) {
    TuiReport report;
    try {
        // Prepare first: a bad option or a mesh that cannot be drawn is reported
        // on a normal screen, before anything about the terminal changes.
        TuiSession session(rMesh, std::move(Options));
        RawTerminal terminal;
        std::string error;
        if (!terminal.Enter(error)) {
            report.mExit = 1;
            report.mError = error;
            return report;
        }
        TerminalIo io(terminal);
        report = session.Run(io);
        terminal.Leave();
    } catch (const std::exception& rErr) {
        report.mExit = 1;
        report.mError = rErr.what();
    }
    return report;
}

}  // namespace meshioplusplus::cli::tui
