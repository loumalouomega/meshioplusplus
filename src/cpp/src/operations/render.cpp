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

// render: the mesh-to-frame half of operations/render.hpp. The tile
// rasterizer is detail/raster.cpp; the text, PNG and snapshot encodings are
// render_text.cpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/render.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/face_color.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/projection.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"
#include "../detail/raster.hpp"

namespace meshioplusplus {
namespace {
constexpr const char* kRndPrefix = "meshio++: render: ";
constexpr double kRndMargin = 0.05;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

void rnd_check_angle(double Value, const char* pName) {
    if (!(Value >= 0.0 && Value <= 180.0))
        throw std::invalid_argument(std::string(kRndPrefix) + pName +
                                    " must lie in [0, 180] degrees");
}

void rnd_validate(const RenderOptions& rOpt) {
    if (rOpt.mWidth <= 0 || rOpt.mHeight <= 0)
        throw std::invalid_argument(std::string(kRndPrefix) + "width and height must be positive");
    if (rOpt.mWidth > 16384 || rOpt.mHeight > 16384)
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "width and height must not exceed 16384");
    if (rOpt.mSupersample != 1 && rOpt.mSupersample != 2 && rOpt.mSupersample != 4)
        throw std::invalid_argument(std::string(kRndPrefix) + "supersample must be 1, 2 or 4");
    if (!(rOpt.mPixelAspect > 0.0) || !std::isfinite(rOpt.mPixelAspect))
        throw std::invalid_argument(std::string(kRndPrefix) + "pixel aspect must be positive");
    if (!(rOpt.mZoom > 0.0) || !std::isfinite(rOpt.mZoom))
        throw std::invalid_argument(std::string(kRndPrefix) + "zoom must be positive");
    if (!(rOpt.mFovDeg > 0.0 && rOpt.mFovDeg < 180.0))
        throw std::invalid_argument(std::string(kRndPrefix) + "fov must lie in (0, 180) degrees");
    if (!(rOpt.mAmbient >= 0.0 && rOpt.mAmbient <= 1.0))
        throw std::invalid_argument(std::string(kRndPrefix) + "ambient must lie in [0, 1]");
    if (!(rOpt.mPointRadius >= 0.0) || !std::isfinite(rOpt.mPointRadius))
        throw std::invalid_argument(std::string(kRndPrefix) + "point radius must be non-negative");
    rnd_check_angle(rOpt.mSplitAngle, "split angle");
    rnd_check_angle(rOpt.mFeatureAngle, "feature angle");
    for (double v : {rOpt.mAzimuth, rOpt.mElevation, rOpt.mRoll, rOpt.mPanX, rOpt.mPanY})
        if (!std::isfinite(v))
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "camera angles and pan must be finite");
    if (rOpt.mVMin.has_value() && rOpt.mVMax.has_value() && *rOpt.mVMin > *rOpt.mVMax)
        throw std::invalid_argument(std::string(kRndPrefix) + "vmin must not exceed vmax");
}

// ---------------------------------------------------------------------------
// Geometry: what is drawn, as rings of point ids into one combined point list
// ---------------------------------------------------------------------------

// One family of drawn primitives that share a colour-resolution context: the
// mesh their local node and cell ids index (the input, or its skin).
struct RndGroup {
    const Mesh* mpDraw = nullptr;
    std::int64_t mPointOffset = 0;  // local point id + offset = combined id
    std::vector<std::int64_t> mStart = {0};
    std::vector<std::int64_t> mNodes;  // local point ids
    std::vector<std::int64_t> mCells;  // local cell index (block-major), per primitive
    std::vector<std::int64_t> mIds;    // input cell index, per primitive (-1 none)
};

struct RndGeometry {
    std::vector<double> mXyz;  // combined points, 3 per point
    std::size_t mNumSource = 0;
    RndGroup mSkinFaces;
    RndGroup mFaces;
    RndGroup mLines;
    RndGroup mPoints;
};

void rnd_append_points(const Mesh& rMesh, std::vector<double>& rXyz) {
    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t base = rXyz.size();
    rXyz.resize(base + 3 * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            rXyz[base + 3 * i + k] = k < dim ? detail::read_double(points, i * dim + k) : 0.0;
}

void rnd_push(RndGroup& rGroup, const std::int64_t* pNodes, std::size_t Count, std::int64_t Cell,
              std::int64_t Id) {
    rGroup.mNodes.insert(rGroup.mNodes.end(), pNodes, pNodes + Count);
    rGroup.mStart.push_back(static_cast<std::int64_t>(rGroup.mNodes.size()));
    rGroup.mCells.push_back(Cell);
    rGroup.mIds.push_back(Id);
}

// The corners a 1-D or 2-D cell type is drawn through; 0 for a type that is not
// drawn this way.
std::size_t rnd_corner_count(const std::string& rType, CellType Type) {
    switch (Type) {
        case CellType::Vertex:
            return 1;
        case CellType::VtkLagrangeCurve:
            return 2;
        case CellType::VtkLagrangeTriangle:
            return 3;
        case CellType::VtkLagrangeQuadrilateral:
            return 4;
        default:
            break;
    }
    if (rType.rfind("line", 0) == 0)
        return 2;
    if (rType.rfind("triangle", 0) == 0)
        return 3;
    if (rType.rfind("quad", 0) == 0)
        return 4;
    return 0;
}

RndGeometry rnd_gather(const Mesh& rMesh, Mesh& rSkin, bool& rHasSkin) {
    RndGeometry g;
    rnd_append_points(rMesh, g.mXyz);
    g.mNumSource = rMesh.NumPoints();
    g.mFaces.mpDraw = &rMesh;
    g.mLines.mpDraw = &rMesh;
    g.mPoints.mpDraw = &rMesh;

    bool has_volume = false;
    bool has_cells = false;
    std::int64_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t block_base = base;
        const std::size_t n = cb.NumCells();
        base += static_cast<std::int64_t>(n);
        if (n > 0)
            has_cells = true;
        const std::string& type = cb.Type();
        const CellType ct = cell_type_from_name(type);
        if (cb.IsPolyhedron() || cell_type_dimension(ct) == 3) {
            has_volume = true;
            continue;
        }
        if (cb.IsRagged()) {
            // A ragged 2-D block is a polygon block: each row is a ring.
            if (cell_type_dimension(ct) != 2 && ct != CellType::Custom) {
                log::warn("render: skipping ragged '{}' block", type);
                continue;
            }
            for (std::size_t r = 0; r < n; ++r) {
                rnd_push(g.mFaces, cb.Row(r), cb.RowSize(r),
                         block_base + static_cast<std::int64_t>(r),
                         block_base + static_cast<std::int64_t>(r));
            }
            continue;
        }
        const std::size_t corners = rnd_corner_count(type, ct);
        const int dim = cell_type_dimension(ct);
        if (corners == 0 || dim < 0) {
            if (n > 0)
                log::warn("render: skipping '{}' cells, which have no drawable form", type);
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t cols = n == 0 ? 0 : conn.Size() / n;
        if (cols < corners)
            continue;
        RndGroup& group = dim == 2 ? g.mFaces : (dim == 1 ? g.mLines : g.mPoints);
        std::int64_t ring[4];
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t k = 0; k < corners; ++k)
                ring[k] = detail::read_int(conn, r * cols + k);
            const std::int64_t cell = block_base + static_cast<std::int64_t>(r);
            rnd_push(group, ring, corners, cell, cell);
        }
    }

    if (has_volume) {
        try {
            rSkin = detail::surface_extract(rMesh, /*forceFaceMode=*/true, /*linearize=*/true,
                                            /*recordParentIds=*/true, "render");
            rHasSkin = true;
        } catch (const std::invalid_argument& rErr) {
            log::warn("render: no boundary skin drawn ({})", rErr.what());
        }
    }
    if (rHasSkin) {
        g.mSkinFaces.mpDraw = &rSkin;
        g.mSkinFaces.mPointOffset = static_cast<std::int64_t>(g.mNumSource);
        rnd_append_points(rSkin, g.mXyz);
        std::vector<std::int64_t> parents;
        if (rSkin.HasCellData("surface:parent_cell"))
            for (std::size_t b = 0; b < rSkin.NumCellBlocks(); ++b) {
                const NDArray& arr = rSkin.CellData("surface:parent_cell", b);
                for (std::size_t i = 0; i < arr.Size(); ++i)
                    parents.push_back(detail::read_int(arr, i));
            }
        std::int64_t skin_base = 0;
        for (const auto cb : rSkin.CellRange()) {
            const std::size_t n = cb.NumCells();
            for (std::size_t r = 0; r < n; ++r) {
                const std::int64_t cell = skin_base + static_cast<std::int64_t>(r);
                const std::int64_t id = static_cast<std::size_t>(cell) < parents.size()
                                            ? parents[static_cast<std::size_t>(cell)]
                                            : -1;
                if (cb.IsRagged()) {
                    rnd_push(g.mSkinFaces, cb.Row(r), cb.RowSize(r), cell, id);
                } else {
                    const NDArray& conn = cb.Conn();
                    const std::size_t cols = conn.Size() / n;
                    std::vector<std::int64_t> ring(cols);
                    for (std::size_t k = 0; k < cols; ++k)
                        ring[k] = detail::read_int(conn, r * cols + k);
                    rnd_push(g.mSkinFaces, ring.data(), cols, cell, id);
                }
            }
            skin_base += static_cast<std::int64_t>(n);
        }
    }

    if (!has_cells) {
        // A cloud: every point is drawn.
        for (std::size_t i = 0; i < rMesh.NumPoints(); ++i) {
            const std::int64_t p = static_cast<std::int64_t>(i);
            rnd_push(g.mPoints, &p, 1, -1, -1);
        }
    }
    return g;
}

// ---------------------------------------------------------------------------
// Colour
// ---------------------------------------------------------------------------

// One value per primitive of a group, by the face_color rules: point data as
// the mean of the corners (summed in node order, divided once), cell data by
// the owning cell. Built on resolve_face_colors so the name lookup, component
// reduction, parent mapping and error messages are the SVG writer's own.
std::vector<double> rnd_group_values(const detail::ColorSpec& rSpec, const Mesh& rSource,
                                     const RndGroup& rGroup, bool IsPoint) {
    const std::size_t count = rGroup.mCells.size();
    std::vector<double> out(count, std::nan(""));
    if (count == 0)
        return out;
    std::vector<detail::ColorFace> faces;
    if (IsPoint) {
        faces.reserve(rGroup.mNodes.size());
        for (std::size_t f = 0; f < count; ++f)
            for (std::int64_t k = rGroup.mStart[f]; k < rGroup.mStart[f + 1]; ++k)
                faces.push_back({{rGroup.mNodes[static_cast<std::size_t>(k)], -1, -1, -1},
                                 1,
                                 rGroup.mCells[f]});
    } else {
        faces.reserve(count);
        for (std::size_t f = 0; f < count; ++f)
            faces.push_back({{-1, -1, -1, -1}, 0, rGroup.mCells[f]});
    }
    const detail::FaceColors fc =
        detail::resolve_face_colors(rSpec, rSource, *rGroup.mpDraw, faces);
    if (!IsPoint) {
        for (std::size_t f = 0; f < count; ++f)
            out[f] = fc.mValues[f];
        return out;
    }
    for (std::size_t f = 0; f < count; ++f) {
        const std::int64_t lo = rGroup.mStart[f];
        const std::int64_t hi = rGroup.mStart[f + 1];
        if (hi == lo)
            continue;
        double sum = 0.0;
        for (std::int64_t k = lo; k < hi; ++k)
            sum += fc.mValues[static_cast<std::size_t>(k)];
        out[f] = sum / static_cast<double>(hi - lo);
    }
    return out;
}

using RndColor = std::array<std::uint8_t, 4>;

struct RndColoring {
    bool mActive = false;
    const std::uint8_t* mpTable = nullptr;
    double mVMin = 0.0;
    double mVMax = 0.0;
    RndColor mNan = {128, 128, 128, 255};

    RndColor Map(double v, const RndColor& rFallback) const {
        if (!mActive)
            return rFallback;
        if (!std::isfinite(v))
            return mNan;
        const detail::Rgb c =
            detail::colormap_lookup(mpTable, detail::color_param(v, mVMin, mVMax));
        return {c.mR, c.mG, c.mB, 255};
    }
};

// ---------------------------------------------------------------------------
// Camera
// ---------------------------------------------------------------------------

double rnd_dot(const double* pA, const std::array<double, 3>& rB) {
    return pA[0] * rB[0] + pA[1] * rB[1] + pA[2] * rB[2];
}

// The Lambert intensity of a unit normal under a unit light direction. A zero
// normal (a degenerate face, a NaN vertex normal) is lit fully rather than
// darkened, since it has no side to turn away.
double rnd_lambert(const double* pNormal, const std::array<double, 3>& rLight,
                   const RenderOptions& rOpt) {
    if (!(std::isfinite(pNormal[0]) && std::isfinite(pNormal[1]) && std::isfinite(pNormal[2])) ||
        (pNormal[0] == 0.0 && pNormal[1] == 0.0 && pNormal[2] == 0.0))
        return 1.0;
    double d = rnd_dot(pNormal, rLight);
    if (rOpt.mTwoSided)
        d = std::fabs(d);
    else if (d < 0.0)
        d = 0.0;
    if (d > 1.0)
        d = 1.0;
    return rOpt.mAmbient + (1.0 - rOpt.mAmbient) * d;
}

// ---------------------------------------------------------------------------
// Overlays, drawn into the final frame after the box filter
// ---------------------------------------------------------------------------

void rnd_plot(Frame& rFrame, std::int64_t X, std::int64_t Y, const RndColor& rColor) {
    if (X < 0 || Y < 0 || X >= rFrame.mWidth || Y >= rFrame.mHeight)
        return;
    const std::size_t idx = static_cast<std::size_t>(Y) * static_cast<std::size_t>(rFrame.mWidth) +
                            static_cast<std::size_t>(X);
    std::copy(rColor.begin(), rColor.end(), rFrame.mRgba.data() + idx * 4);
    rFrame.mCellIds[idx] = -1;
}

void rnd_line(Frame& rFrame, double X0, double Y0, double X1, double Y1, const RndColor& rColor) {
    const double span = std::max(std::fabs(X1 - X0), std::fabs(Y1 - Y0));
    const std::int64_t n = std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(span)));
    for (std::int64_t i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n);
        rnd_plot(rFrame, static_cast<std::int64_t>(std::floor(X0 + t * (X1 - X0))),
                 static_cast<std::int64_t>(std::floor(Y0 + t * (Y1 - Y0))), rColor);
    }
}

std::string rnd_num(double v) {
    char buf[64];
    detail::snprintf_c(buf, sizeof(buf), "%.6g", v);
    return buf;
}

// The largest 1, 2 or 5 times a power of ten not above Target.
double rnd_nice_length(double Target) {
    const double base = std::pow(10.0, std::floor(std::log10(Target)));
    for (double m : {5.0, 2.0, 1.0})
        if (m * base <= Target)
            return m * base;
    return base;
}

}  // namespace

namespace detail {

// The camera direction for a named view: the camera sits on that side.
void render_view_angles(const RenderOptions& rOpt, double& rAzimuth, double& rElevation) {
    rAzimuth = rOpt.mAzimuth;
    rElevation = rOpt.mElevation;
    const std::string& v = rOpt.mView;
    if (v.empty())
        return;
    if (v == "iso") {
        rAzimuth = 45.0;
        rElevation = 35.264389682754654;
    } else if (v == "+x") {
        rAzimuth = 0.0;
        rElevation = 0.0;
    } else if (v == "-x") {
        rAzimuth = 180.0;
        rElevation = 0.0;
    } else if (v == "+y") {
        rAzimuth = 90.0;
        rElevation = 0.0;
    } else if (v == "-y") {
        rAzimuth = -90.0;
        rElevation = 0.0;
    } else if (v == "+z") {
        rAzimuth = -90.0;
        rElevation = 90.0;
    } else if (v == "-z") {
        rAzimuth = -90.0;
        rElevation = -90.0;
    } else {
        throw std::invalid_argument(std::string(kRndPrefix) + "unknown view '" + v +
                                    "' (expected iso, +x, -x, +y, -y, +z or -z)");
    }
}

}  // namespace detail

Frame render(const Mesh& rMesh, const RenderOptions& rOpt) {
    rnd_validate(rOpt);
    double azimuth = 0.0;
    double elevation = 0.0;
    detail::render_view_angles(rOpt, azimuth, elevation);
    const detail::CameraBasis cam = detail::camera_basis(azimuth, elevation, rOpt.mRoll);

    Mesh skin;
    bool has_skin = false;
    RndGeometry g = rnd_gather(rMesh, skin, has_skin);

    // The drawn faces as rings into the combined points: the skin's first, then
    // the input's own 2-D cells.
    std::vector<std::int64_t> f_start = {0};
    std::vector<std::int64_t> f_nodes;
    std::vector<std::int64_t> f_ids;
    for (const RndGroup* p_group : {&g.mSkinFaces, &g.mFaces}) {
        const RndGroup& group = *p_group;
        for (std::size_t f = 0; f < group.mIds.size(); ++f) {
            for (std::int64_t k = group.mStart[f]; k < group.mStart[f + 1]; ++k)
                f_nodes.push_back(group.mNodes[static_cast<std::size_t>(k)] + group.mPointOffset);
            f_start.push_back(static_cast<std::int64_t>(f_nodes.size()));
            f_ids.push_back(group.mIds[f]);
        }
    }
    const std::size_t num_faces = f_ids.size();

    // Values and the mapped range.
    RndColoring coloring;
    coloring.mNan = rOpt.mNanColor;
    std::vector<double> f_values(num_faces, 0.0);
    std::vector<double> l_values(g.mLines.mIds.size(), 0.0);
    std::vector<double> p_values(g.mPoints.mIds.size(), 0.0);
    if (!rOpt.mColorBy.empty()) {
        detail::ColorSpec spec;
        spec.mColorBy = rOpt.mColorBy;
        spec.mComponent = rOpt.mComponent;
        spec.mCmap = rOpt.mCmap;
        coloring.mpTable = detail::colormap_table(rOpt.mCmap);
        const bool is_point = rMesh.HasPointData(rOpt.mColorBy);
        if (!is_point && !rMesh.HasCellData(rOpt.mColorBy))
            (void)detail::resolve_face_colors(spec, rMesh, rMesh, {});  // throws the shared message
        f_values = rnd_group_values(spec, rMesh, g.mSkinFaces, is_point);
        const std::vector<double> own = rnd_group_values(spec, rMesh, g.mFaces, is_point);
        f_values.insert(f_values.end(), own.begin(), own.end());
        l_values = rnd_group_values(spec, rMesh, g.mLines, is_point);
        p_values = rnd_group_values(spec, rMesh, g.mPoints, is_point);
        double lo = 0.0;
        double hi = 0.0;
        bool seen = false;
        for (const std::vector<double>* p_values_list : {&f_values, &l_values, &p_values}) {
            for (double v : *p_values_list) {
                if (!std::isfinite(v))
                    continue;
                if (!seen) {
                    lo = v;
                    hi = v;
                    seen = true;
                } else {
                    lo = std::min(lo, v);
                    hi = std::max(hi, v);
                }
            }
        }
        coloring.mVMin = rOpt.mVMin.value_or(lo);
        coloring.mVMax = rOpt.mVMax.value_or(hi);
        if (coloring.mVMin > coloring.mVMax)
            throw std::invalid_argument(std::string(kRndPrefix) + "vmin must not exceed vmax");
        coloring.mActive = true;
    }

    // Face normals (Newell's, from the shared crease kernel).
    std::vector<double> f_normal(3 * num_faces, 0.0);
    parallel_for(num_faces, [&](std::size_t f) {
        detail::ring_unit_normal(g.mXyz.data(), f_nodes.data() + f_start[f],
                                 static_cast<std::size_t>(f_start[f + 1] - f_start[f]),
                                 &f_normal[3 * f]);
    });

    // Smooth shading draws the faces through compute_normals' split points:
    // the rings change, the original points keep their indices.
    std::vector<double> xyz = std::move(g.mXyz);
    std::vector<std::int64_t> s_nodes = f_nodes;
    std::vector<double> v_normal;
    bool smooth = rOpt.mShading == RenderShading::Smooth && num_faces > 0;
    if (smooth) {
        try {
            Mesh poly;
            const std::size_t nv = xyz.size() / 3;
            NDArray pts = NDArray::Uninit(DType::Float64, {nv, 3});
            std::copy(xyz.begin(), xyz.end(), pts.As<double>());
            poly.AssignPoints(std::move(pts));
            poly.AddPolygonBlock("polygon", std::vector<std::int64_t>(f_nodes),
                                 std::vector<std::int64_t>(f_start));
            NormalsOptions normals_options;
            normals_options.mSplit = true;
            normals_options.mSplitAngle = rOpt.mSplitAngle;
            const NormalsResult result = compute_normals(poly, normals_options);
            const Mesh& m = result.mMesh;
            std::vector<double> split_xyz;
            rnd_append_points(m, split_xyz);
            const auto block = m.Cells(0);
            for (std::size_t f = 0; f < num_faces; ++f) {
                const std::int64_t* row = block.Row(f);
                for (std::int64_t k = f_start[f]; k < f_start[f + 1]; ++k)
                    s_nodes[static_cast<std::size_t>(k)] = row[k - f_start[f]];
            }
            const NDArray& normals = m.PointData(kNormalsName);
            v_normal.resize(normals.Size());
            for (std::size_t i = 0; i < normals.Size(); ++i)
                v_normal[i] = detail::read_double(normals, i);
            xyz = std::move(split_xyz);
        } catch (const std::exception& rErr) {
            log::warn("render: smooth shading fell back to flat ({})", rErr.what());
            smooth = false;
            s_nodes = f_nodes;
        }
    }

    // Project every vertex a primitive uses.
    const std::size_t nv = xyz.size() / 3;
    std::vector<std::uint8_t> used(nv, 0);
    for (std::int64_t v : s_nodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : g.mLines.mNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : g.mPoints.mNodes)
        used[static_cast<std::size_t>(v)] = 1;

    const bool perspective = rOpt.mProjection == RenderProjection::Perspective;
    std::array<double, 3> eye = {0.0, 0.0, 0.0};
    const double focal = 1.0 / std::tan(rOpt.mFovDeg * 3.141592653589793 / 360.0);
    if (perspective) {
        // Fit the bounding sphere of the drawn vertices into the field of view.
        std::array<double, 3> lo = {0.0, 0.0, 0.0};
        std::array<double, 3> hi = {0.0, 0.0, 0.0};
        bool seen = false;
        for (std::size_t i = 0; i < nv; ++i) {
            if (!used[i])
                continue;
            for (int k = 0; k < 3; ++k) {
                const double c = xyz[3 * i + static_cast<std::size_t>(k)];
                lo[k] = seen ? std::min(lo[k], c) : c;
                hi[k] = seen ? std::max(hi[k], c) : c;
            }
            seen = true;
        }
        std::array<double, 3> centre;
        for (int k = 0; k < 3; ++k)
            centre[k] = 0.5 * (lo[k] + hi[k]);
        double radius = 0.0;
        for (std::size_t i = 0; i < nv; ++i) {
            if (!used[i])
                continue;
            double d2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double d = xyz[3 * i + static_cast<std::size_t>(k)] - centre[k];
                d2 += d * d;
            }
            radius = std::max(radius, std::sqrt(d2));
        }
        if (!(radius > 0.0))
            radius = 1.0;
        const double distance = radius / std::sin(rOpt.mFovDeg * 3.141592653589793 / 360.0);
        for (int k = 0; k < 3; ++k)
            eye[k] = centre[k] + cam.mW[static_cast<std::size_t>(k)] * distance;
    }

    std::vector<double> sx(nv, 0.0);
    std::vector<double> sy(nv, 0.0);
    std::vector<double> sz(nv, 0.0);
    parallel_for(nv, [&](std::size_t i) {
        const double* p = &xyz[3 * i];
        if (perspective) {
            const double rel[3] = {p[0] - eye[0], p[1] - eye[1], p[2] - eye[2]};
            const double zc = -rnd_dot(rel, cam.mW);
            sx[i] = focal * rnd_dot(rel, cam.mU) / zc;
            sy[i] = focal * rnd_dot(rel, cam.mV) / zc;
            sz[i] = 1.0 / zc;
        } else {
            sx[i] = rnd_dot(p, cam.mU);
            sy[i] = rnd_dot(p, cam.mV);
            sz[i] = rnd_dot(p, cam.mW);
        }
    });

    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0, min_z = 0.0, max_z = 0.0;
    bool any = false;
    for (std::size_t i = 0; i < nv; ++i) {
        if (!used[i] || !std::isfinite(sx[i]) || !std::isfinite(sy[i]) || !std::isfinite(sz[i]))
            continue;
        if (!any) {
            min_x = max_x = sx[i];
            min_y = max_y = sy[i];
            min_z = max_z = sz[i];
            any = true;
        } else {
            min_x = std::min(min_x, sx[i]);
            max_x = std::max(max_x, sx[i]);
            min_y = std::min(min_y, sy[i]);
            max_y = std::max(max_y, sy[i]);
            min_z = std::min(min_z, sz[i]);
            max_z = std::max(max_z, sz[i]);
        }
    }

    // Fit the drawn geometry to the frame (leaving room for a colour bar).
    const int ss = rOpt.mSupersample;
    const double width = static_cast<double>(rOpt.mWidth) * ss;
    const double height = static_cast<double>(rOpt.mHeight) * ss;
    const double aspect = rOpt.mPixelAspect;
    const bool colorbar = rOpt.mColorbar && coloring.mActive;
    const double usable_w = width * (colorbar ? 0.85 : 1.0);
    const double avail_x = usable_w * (1.0 - 2.0 * kRndMargin);
    const double avail_y = height * aspect * (1.0 - 2.0 * kRndMargin);
    const double extent_x = max_x - min_x;
    const double extent_y = max_y - min_y;
    double scale = std::numeric_limits<double>::infinity();
    if (extent_x > 0.0)
        scale = std::min(scale, avail_x / extent_x);
    if (extent_y > 0.0)
        scale = std::min(scale, avail_y / extent_y);
    if (!std::isfinite(scale))
        scale = std::min(avail_x, avail_y);
    scale *= rOpt.mZoom;
    const double mid_x = 0.5 * (min_x + max_x);
    const double mid_y = 0.5 * (min_y + max_y);
    const double depth_range = max_z - min_z;

    detail::RasterScene scene;
    scene.mVertices.resize(nv);
    parallel_for(nv, [&](std::size_t i) {
        detail::RasterVertex& v = scene.mVertices[i];
        v.mX = 0.5 * usable_w + (sx[i] - mid_x) * scale + rOpt.mPanX * width;
        v.mY = 0.5 * height - (sy[i] - mid_y) * scale / aspect - rOpt.mPanY * height;
        v.mDepth = depth_range > 0.0 ? (sz[i] - min_z) / depth_range : 0.5;
    });
    scene.mLineWidth = ss;
    scene.mPointRadius = rOpt.mPointRadius * ss;
    // A line samples depth along its own path, a face at pixel centres: on a
    // face seen at an angle the two differ by the depth slope over a pixel or
    // so. The normalized depth spans roughly the frame, so three pixels' worth
    // of it (plus a floor) keeps an edge on top of its own face while anything
    // a real distance behind still loses.
    scene.mDepthBias = 1.0e-3 + 3.0 / std::max(width, height);

    // The light, from camera space into world space.
    std::array<double, 3> light;
    for (int k = 0; k < 3; ++k)
        light[k] = rOpt.mLightDir[0] * cam.mU[k] + rOpt.mLightDir[1] * cam.mV[k] +
                   rOpt.mLightDir[2] * cam.mW[k];
    const double light_len =
        std::sqrt(light[0] * light[0] + light[1] * light[1] + light[2] * light[2]);
    if (light_len > 0.0)
        for (double& c : light)
            c /= light_len;
    else
        light = cam.mW;

    // Triangles: a triangle as is, a quad split on its shorter diagonal, a
    // polygon fanned from its first corner.
    std::size_t num_tris = 0;
    for (std::size_t f = 0; f < num_faces; ++f) {
        const std::int64_t n = f_start[f + 1] - f_start[f];
        if (n >= 3)
            num_tris += static_cast<std::size_t>(n - 2);
    }
    scene.mTris.reserve(num_tris);
    auto dist2 = [&](std::int64_t a, std::int64_t b) {
        double d2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const double d =
                xyz[3 * static_cast<std::size_t>(a) + k] - xyz[3 * static_cast<std::size_t>(b) + k];
            d2 += d * d;
        }
        return d2;
    };
    for (std::size_t f = 0; f < num_faces; ++f) {
        const std::int64_t* ring = s_nodes.data() + f_start[f];
        const std::int64_t n = f_start[f + 1] - f_start[f];
        if (n < 3)
            continue;
        const RndColor color = coloring.Map(f_values[f], rOpt.mFillColor);
        const double face_i =
            rOpt.mShading == RenderShading::None ? 1.0 : rnd_lambert(&f_normal[3 * f], light, rOpt);
        auto corner_i = [&](std::int64_t v) {
            if (!smooth)
                return face_i;
            const double* nrm = &v_normal[3 * static_cast<std::size_t>(v)];
            if (!(std::isfinite(nrm[0]) && std::isfinite(nrm[1]) && std::isfinite(nrm[2])))
                return face_i;
            return rnd_lambert(nrm, light, rOpt);
        };
        auto emit = [&](std::int64_t a, std::int64_t b, std::int64_t c) {
            detail::RasterTri t;
            t.mV[0] = ring[a];
            t.mV[1] = ring[b];
            t.mV[2] = ring[c];
            t.mIntensity[0] = corner_i(ring[a]);
            t.mIntensity[1] = corner_i(ring[b]);
            t.mIntensity[2] = corner_i(ring[c]);
            t.mColor = color;
            t.mId = f_ids[f];
            scene.mTris.push_back(t);
        };
        if (n == 4) {
            if (dist2(ring[0], ring[2]) <= dist2(ring[1], ring[3])) {
                emit(0, 1, 2);
                emit(0, 2, 3);
            } else {
                emit(1, 2, 3);
                emit(1, 3, 0);
            }
        } else {
            for (std::int64_t k = 1; k + 1 < n; ++k)
                emit(0, k, k + 1);
        }
    }

    // Lines: the input's line cells, then the edge overlay.
    for (std::size_t i = 0; i < g.mLines.mIds.size(); ++i) {
        detail::RasterLine line;
        line.mA = g.mLines.mNodes[static_cast<std::size_t>(g.mLines.mStart[i])];
        line.mB = g.mLines.mNodes[static_cast<std::size_t>(g.mLines.mStart[i]) + 1];
        line.mColor = coloring.Map(l_values[i], rOpt.mLineColor);
        line.mId = g.mLines.mIds[i];
        scene.mLines.push_back(line);
    }
    if (rOpt.mEdges == RenderEdges::All) {
        std::vector<std::pair<std::int64_t, std::int64_t>> edges;
        edges.reserve(f_nodes.size());
        for (std::size_t f = 0; f < num_faces; ++f) {
            const std::int64_t lo = f_start[f];
            const std::int64_t n = f_start[f + 1] - lo;
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t a = f_nodes[static_cast<std::size_t>(lo + k)];
                const std::int64_t b = f_nodes[static_cast<std::size_t>(lo + (k + 1) % n)];
                if (a != b)
                    edges.emplace_back(std::min(a, b), std::max(a, b));
            }
        }
        parallel_sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        for (const auto& e : edges)
            scene.mLines.push_back({e.first, e.second, rOpt.mEdgeColor, -1});
    } else if (rOpt.mEdges == RenderEdges::Feature && num_faces > 0) {
        const std::vector<detail::CreaseEdge> edges =
            detail::crease_edges(f_start, f_nodes, f_normal, rOpt.mFeatureAngle);
        for (const detail::CreaseEdge& e : edges)
            if (e.IsCrease() || e.IsBoundary())
                scene.mLines.push_back({e.mLo, e.mHi, rOpt.mEdgeColor, -1});
    }
    for (std::size_t i = 0; i < g.mPoints.mIds.size(); ++i) {
        detail::RasterPoint point;
        point.mV = g.mPoints.mNodes[static_cast<std::size_t>(g.mPoints.mStart[i])];
        point.mColor = coloring.Map(p_values[i], rOpt.mLineColor);
        point.mId = g.mPoints.mIds[i];
        scene.mPoints.push_back(point);
    }

    const detail::RasterTarget target =
        detail::rasterize(scene, rOpt.mWidth * ss, rOpt.mHeight * ss, rOpt.mBackground);
    Frame frame;
    frame.mWidth = rOpt.mWidth;
    frame.mHeight = rOpt.mHeight;
    detail::raster_downsample(target, ss, frame.mRgba, frame.mCellIds);

    // Overlays and notes, at the output resolution.
    const double min_side = static_cast<double>(std::min(rOpt.mWidth, rOpt.mHeight));
    const std::int64_t margin =
        std::max<std::int64_t>(1, static_cast<std::int64_t>(min_side * 0.03));
    std::int64_t bar_w = 0;
    if (coloring.mActive) {
        frame.mColored = true;
        frame.mVMin = coloring.mVMin;
        frame.mVMax = coloring.mVMax;
        std::string label = rOpt.mColorBy;
        if (rOpt.mComponent.has_value())
            label += "[" + std::to_string(*rOpt.mComponent) + "]";
        frame.mNotes.push_back(label + ": " + rnd_num(coloring.mVMin) + " .. " +
                               rnd_num(coloring.mVMax) + " (" + rOpt.mCmap + ")");
    }
    if (colorbar) {
        bar_w = std::max<std::int64_t>(2, rOpt.mWidth / 40);
        const std::int64_t x1 = rOpt.mWidth - 1 - margin;
        const std::int64_t x0 = x1 - bar_w + 1;
        const std::int64_t y0 = static_cast<std::int64_t>(rOpt.mHeight * 0.2);
        const std::int64_t y1 =
            std::max<std::int64_t>(y0 + 1, static_cast<std::int64_t>(rOpt.mHeight * 0.8));
        for (std::int64_t y = y0; y <= y1; ++y) {
            const double t = static_cast<double>(y1 - y) / static_cast<double>(y1 - y0);
            const detail::Rgb c = detail::colormap_lookup(coloring.mpTable, t);
            for (std::int64_t x = x0; x <= x1; ++x)
                rnd_plot(frame, x, y, {c.mR, c.mG, c.mB, 255});
        }
    }
    if (rOpt.mAxes) {
        const double len = std::max(4.0, min_side * 0.12);
        const double ox = static_cast<double>(margin) + len;
        const double oy = static_cast<double>(rOpt.mHeight - 1 - margin) - len;
        const RndColor axis_colors[3] = {
            {230, 60, 60, 255}, {60, 180, 75, 255}, {70, 110, 240, 255}};
        int order[3] = {0, 1, 2};
        // Farthest axis first, so the nearer ones are drawn over it.
        std::stable_sort(order, order + 3, [&](int a, int b) { return cam.mW[a] < cam.mW[b]; });
        for (int k : order)
            rnd_line(frame, ox, oy, ox + cam.mU[k] * len, oy - cam.mV[k] * len / aspect,
                     axis_colors[k]);
    }
    if (rOpt.mScaleBar) {
        if (perspective) {
            frame.mNotes.push_back("scale bar: not drawn with a perspective camera");
        } else if (any) {
            const double px_per_unit = scale / ss;
            const double length = rnd_nice_length(0.2 * rOpt.mWidth / px_per_unit);
            const double length_px = length * px_per_unit;
            const double x1 =
                static_cast<double>(rOpt.mWidth - 1 - margin - (bar_w > 0 ? bar_w + margin : 0));
            const double x0 = x1 - length_px;
            const double y = static_cast<double>(rOpt.mHeight - 1 - margin);
            const RndColor grey = {128, 128, 128, 255};
            rnd_line(frame, x0, y, x1, y, grey);
            rnd_line(frame, x0, y - 2, x0, y, grey);
            rnd_line(frame, x1, y - 2, x1, y, grey);
            frame.mNotes.push_back("scale bar: " + rnd_num(length));
        }
    }
    return frame;
}

}  // namespace meshioplusplus
