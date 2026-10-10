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
#include <memory>
#include <optional>
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
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/tensor_invariants.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"
#include "../detail/raster.hpp"
#include "../detail/render_field.hpp"
#include "../detail/streamlines.hpp"

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
    if (rOpt.mClipLow.has_value() || rOpt.mClipHigh.has_value()) {
        const double lo = rOpt.mClipLow.value_or(0.0);
        const double hi = rOpt.mClipHigh.value_or(100.0);
        if (!(lo >= 0.0 && hi <= 100.0 && lo < hi))
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "clip percentiles must satisfy 0 <= low < high <= 100");
    }
    if (!(rOpt.mScaleThreshold > 0.0) || !std::isfinite(rOpt.mScaleThreshold))
        throw std::invalid_argument(std::string(kRndPrefix) + "scale threshold must be positive");
    if (rOpt.mIsolines < 0 || rOpt.mIsolines > 1000)
        throw std::invalid_argument(std::string(kRndPrefix) + "isolines must lie in [0, 1000]");
    for (double level : rOpt.mIsoLevels)
        if (!std::isfinite(level))
            throw std::invalid_argument(std::string(kRndPrefix) + "isoline levels must be finite");
    if (rOpt.mVectorCount < 1 || rOpt.mVectorCount > 100000)
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "vector count must lie in [1, 100000]");
    if (rOpt.mGridFloor && rOpt.mTheme == RenderTheme::None)
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "the grid floor is part of a theme (set theme)");
    if (rOpt.mStreamSeeds < 1 || rOpt.mStreamSeeds > 10000)
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "stream seeds must lie in [1, 10000]");
    if (!(rOpt.mStreamLength > 0.0) || !(rOpt.mStreamLength <= 10.0))
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "stream length must lie in (0, 10] diagonals");
    if (!(rOpt.mVectorLength >= 0.0) || !std::isfinite(rOpt.mVectorLength) ||
        !std::isfinite(rOpt.mWarpScale))
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "vector length must be non-negative and the warp scale finite");
    if (!rOpt.mExpr.empty() && !rOpt.mColorBy.empty())
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "give either color_by or expr, not both");
    if (!rOpt.mReduce.empty() && rOpt.mColorBy.empty())
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "reduce needs color_by (a tensor array)");
    if (rOpt.mColorRegions && (!rOpt.mColorBy.empty() || !rOpt.mExpr.empty()))
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "color_regions excludes color_by and expr");
    if (rOpt.mDiagnostic != RenderDiagnostic::None &&
        (!rOpt.mColorBy.empty() || !rOpt.mExpr.empty() || rOpt.mColorRegions))
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "a diagnostic view excludes color_by, expr and color_regions");
    if (rOpt.mDiagnostic == RenderDiagnostic::Quality && rOpt.mQualityMetric.empty())
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "the quality diagnostic needs a metric (quality_metric)");
    if (rOpt.mCutaways.size() > 2)
        throw std::invalid_argument(std::string(kRndPrefix) + "at most two cutaway planes");
    for (const RenderCutaway& cut : rOpt.mCutaways) {
        double len2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            if (!std::isfinite(cut.mPoint[k]) || !std::isfinite(cut.mNormal[k]))
                throw std::invalid_argument(std::string(kRndPrefix) +
                                            "a cutaway plane must be finite");
            len2 += cut.mNormal[k] * cut.mNormal[k];
        }
        if (!(len2 > 0.0))
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "a cutaway plane needs a non-zero normal");
    }
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
    RenderScale mScale = RenderScale::Linear;
    double mThreshold = 1.0;
    /// Categories: a value is an index into this palette (cycling).
    std::vector<RndColor> mPalette;
    /// A diagnostic that colours faces only: lines and points keep their own.
    bool mFacesOnly = false;

    RndColor Map(double v, const RndColor& rFallback) const {
        if (!mActive)
            return rFallback;
        if (!mPalette.empty()) {
            if (!std::isfinite(v) || v < 0.0)
                return mNan;
            return mPalette[static_cast<std::size_t>(v) % mPalette.size()];
        }
        const double s = detail::render_scale_forward(mScale, mThreshold, v);
        if (!std::isfinite(s))
            return mNan;
        const double lo = detail::render_scale_forward(mScale, mThreshold, mVMin);
        const double hi = detail::render_scale_forward(mScale, mThreshold, mVMax);
        const detail::Rgb c = detail::colormap_lookup(mpTable, detail::color_param(s, lo, hi));
        return {c.mR, c.mG, c.mB, 255};
    }

    /// The colour of a line or point, which a faces-only diagnostic leaves alone.
    RndColor MapOther(double v, const RndColor& rFallback) const {
        return mFacesOnly ? rFallback : Map(v, rFallback);
    }
};

std::string rnd_hex(const RndColor& rColor) {
    static const char digits[] = "0123456789abcdef";
    std::string out = "#";
    for (int k = 0; k < 3; ++k) {
        out.push_back(digits[(rColor[static_cast<std::size_t>(k)] >> 4) & 15]);
        out.push_back(digits[rColor[static_cast<std::size_t>(k)] & 15]);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Derived arrays and per-vertex data
// ---------------------------------------------------------------------------

std::string rnd_available(const Mesh& rMesh) {
    std::string names;
    for (const std::string& n : rMesh.PointDataNames())
        names += (names.empty() ? "" : ", ") + n;
    return names.empty() ? "none" : names;
}

// The array the colouring reads, when it is not simply `mColorBy`: a
// `data_calc` expression, a tensor invariant, or a `quality` metric. They all
// become an ordinary array on a copy of the input, so everything after this
// point -- the skin, the parent cell mapping, the shared resolver -- is the
// path an array of the input takes.
void rnd_prepare_source(const RenderOptions& rOpt, const Mesh& rMesh, Mesh& rWork,
                        const Mesh*& rpSrc, std::string& rColorBy, std::string& rLabel,
                        std::optional<int>& rComponent) {
    rpSrc = &rMesh;
    rColorBy = rOpt.mColorBy;
    rLabel = rColorBy;
    rComponent = rOpt.mComponent;
    if (!rOpt.mExpr.empty()) {
        DataCalcOptions options;
        options.output = "render:expr";
        options.overwrite = true;
        try {
            options.location = DataLocation::Point;
            rWork = data_calc(rMesh, rOpt.mExpr, options);
        } catch (const std::invalid_argument& rAsPoint) {
            try {
                options.location = DataLocation::Cell;
                rWork = data_calc(rMesh, rOpt.mExpr, options);
            } catch (const std::invalid_argument& rAsCell) {
                throw std::invalid_argument(std::string(kRndPrefix) + "expression '" + rOpt.mExpr +
                                            "' failed as point data (" + rAsPoint.what() +
                                            ") and as cell data (" + rAsCell.what() + ")");
            }
        }
        rpSrc = &rWork;
        rColorBy = "render:expr";
        rLabel = rOpt.mExpr;
        return;
    }
    if (!rOpt.mReduce.empty()) {
        if (rOpt.mReduce != "mises" && rOpt.mReduce != "hydrostatic" && rOpt.mReduce != "principal")
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "reduce must be 'mises', "
                                        "'hydrostatic' or 'principal', not '" +
                                        rOpt.mReduce + "'");
        TensorInvariantsOptions options;
        if (rMesh.HasPointData(rColorBy))
            options.location = DataLocation::Point;
        else if (rMesh.HasCellData(rColorBy))
            options.location = DataLocation::Cell;
        else
            throw std::invalid_argument(std::string(kRndPrefix) + "no array named '" + rColorBy +
                                        "' to reduce (point data: " + rnd_available(rMesh) + ")");
        options.names = {rColorBy};
        options.prefix = "render:";
        options.outputs = rOpt.mReduce == "mises"       ? TensorInvariant::Mises
                          : rOpt.mReduce == "principal" ? TensorInvariant::Principal
                                                        : TensorInvariant::Hydrostatic;
        rWork = tensor_invariants(rMesh, options);
        rpSrc = &rWork;
        rColorBy = "render:" + rColorBy + "_" + rOpt.mReduce;
        rLabel = rOpt.mReduce + "(" + rOpt.mColorBy + ")";
        if (rOpt.mReduce == "principal" && !rComponent.has_value())
            rComponent = 2;
        else if (rOpt.mReduce != "principal")
            rComponent.reset();
        return;
    }
    if (rOpt.mDiagnostic == RenderDiagnostic::Quality ||
        rOpt.mDiagnostic == RenderDiagnostic::Inverted ||
        rOpt.mDiagnostic == RenderDiagnostic::Degenerate) {
        rWork = attach_quality(rMesh);
        rpSrc = &rWork;
        rColorBy = rOpt.mDiagnostic == RenderDiagnostic::Quality
                       ? "quality:" + rOpt.mQualityMetric
                       : (rOpt.mDiagnostic == RenderDiagnostic::Inverted ? "quality:inverted"
                                                                         : "quality:degenerate");
        rLabel = rColorBy;
    }
}

// One value per combined vertex of a point array, read as `Ncomp` components
// per vertex (missing ones zero), with the skin's points after the source's.
std::vector<double> rnd_vertex_array(const Mesh& rSrc, const Mesh* pSkin, std::size_t NumSource,
                                     const std::string& rName, std::size_t Ncomp,
                                     const char* pWhat) {
    if (!rSrc.HasPointData(rName))
        throw std::invalid_argument(std::string(kRndPrefix) + pWhat + " array '" + rName +
                                    "' is not point data (available: " + rnd_available(rSrc) + ")");
    const std::size_t total = NumSource + (pSkin != nullptr ? pSkin->NumPoints() : 0);
    std::vector<double> out(Ncomp * total, 0.0);
    auto read = [&](const Mesh& rMesh, std::size_t Offset) {
        const NDArray& arr = rMesh.PointData(rName);
        const std::size_t n = rMesh.NumPoints();
        const std::size_t have = n == 0 ? 0 : arr.Size() / n;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < Ncomp && k < have; ++k)
                out[(Offset + i) * Ncomp + k] = detail::read_double(arr, i * have + k);
    };
    read(rSrc, 0);
    if (pSkin != nullptr && pSkin->HasPointData(rName))
        read(*pSkin, NumSource);
    return out;
}

// The scalar of a point array at every combined vertex, by the rules the face
// colouring uses (a component, or the magnitude).
std::vector<double> rnd_vertex_scalar(const detail::ColorSpec& rSpec, const Mesh& rSrc,
                                      const Mesh* pSkin, std::size_t NumSource) {
    std::vector<detail::ColorFace> faces;
    faces.reserve(rSrc.NumPoints());
    for (std::size_t i = 0; i < rSrc.NumPoints(); ++i)
        faces.push_back({{static_cast<std::int64_t>(i), -1, -1, -1}, 1, -1});
    std::vector<double> out = detail::resolve_face_colors(rSpec, rSrc, rSrc, faces).mValues;
    if (pSkin != nullptr) {
        faces.clear();
        for (std::size_t i = 0; i < pSkin->NumPoints(); ++i)
            faces.push_back({{static_cast<std::int64_t>(i), -1, -1, -1}, 1, -1});
        const std::vector<double> more =
            detail::resolve_face_colors(rSpec, rSrc, *pSkin, faces).mValues;
        out.insert(out.end(), more.begin(), more.end());
    }
    out.resize(NumSource + (pSkin != nullptr ? pSkin->NumPoints() : 0), std::nan(""));
    return out;
}

// Triangles of a face, as positions in its ring: a triangle itself, a quad
// split on its shorter diagonal, a polygon fanned from its first corner.
void rnd_face_triangles(const std::vector<double>& rXyz, const std::int64_t* pRing, std::int64_t N,
                        std::vector<std::array<int, 3>>& rOut) {
    rOut.clear();
    if (N < 3)
        return;
    auto dist2 = [&](std::int64_t a, std::int64_t b) {
        double d2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const double d = rXyz[3 * static_cast<std::size_t>(a) + k] -
                             rXyz[3 * static_cast<std::size_t>(b) + k];
            d2 += d * d;
        }
        return d2;
    };
    if (N == 4) {
        if (dist2(pRing[0], pRing[2]) <= dist2(pRing[1], pRing[3])) {
            rOut.push_back({0, 1, 2});
            rOut.push_back({0, 2, 3});
        } else {
            rOut.push_back({1, 2, 3});
            rOut.push_back({1, 3, 0});
        }
    } else {
        for (int k = 1; k + 1 < N; ++k)
            rOut.push_back({0, k, k + 1});
    }
}

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


// ---------------------------------------------------------------------------
// Theme: a banded background, a grid floor and three post-processes
// ---------------------------------------------------------------------------

int rnd_lerp_i(int A, int B, int Num, int Den) { return A + (B - A) * Num / Den; }

// One colour of a banded gradient through `Stops`: the band `Band` of `Bands`.
RndColor rnd_band_color(const std::vector<std::array<int, 3>>& rStops, int Band, int Bands) {
    const int segments = static_cast<int>(rStops.size()) - 1;
    const int scaled = Bands > 1 ? Band * segments * 1000 / (Bands - 1) : 0;
    int seg = std::min(segments - 1, scaled / 1000);
    const int frac = scaled - seg * 1000;
    const std::array<int, 3>& a = rStops[static_cast<std::size_t>(seg)];
    const std::array<int, 3>& b = rStops[static_cast<std::size_t>(seg) + 1];
    return {static_cast<std::uint8_t>(rnd_lerp_i(a[0], b[0], frac, 1000)),
            static_cast<std::uint8_t>(rnd_lerp_i(a[1], b[1], frac, 1000)),
            static_cast<std::uint8_t>(rnd_lerp_i(a[2], b[2], frac, 1000)), 255};
}

// Put the theme behind the model: a sunset in eight bands down to the horizon
// at 60% of the height, a dark floor in four below it, the grid on the floor
// when asked for; then lay the frame over it by its own alpha.
void rnd_theme_background(Frame& rFrame, const RenderOptions& rOpt) {
    const std::int64_t w = rFrame.mWidth;
    const std::int64_t h = rFrame.mHeight;
    if (w < 1 || h < 1)
        return;
    const std::int64_t horizon = h * 6 / 10;
    static const std::vector<std::array<int, 3>> sky = {
        {20, 8, 60}, {120, 20, 140}, {255, 45, 150}, {255, 150, 40}};
    static const std::vector<std::array<int, 3>> floor_stops = {{12, 4, 36}, {26, 8, 60}};
    Frame bg;
    bg.mWidth = static_cast<int>(w);
    bg.mHeight = static_cast<int>(h);
    bg.mRgba.resize(static_cast<std::size_t>(w * h) * 4);
    bg.mCellIds.assign(static_cast<std::size_t>(w * h), -1);
    for (std::int64_t y = 0; y < h; ++y) {
        const RndColor c = y < horizon
                               ? rnd_band_color(sky, static_cast<int>(y * 8 / std::max<std::int64_t>(1, horizon)), 8)
                               : rnd_band_color(floor_stops,
                                                static_cast<int>((y - horizon) * 4 /
                                                                 std::max<std::int64_t>(1, h - horizon)),
                                                4);
        for (std::int64_t x = 0; x < w; ++x)
            std::copy(c.begin(), c.end(), bg.mRgba.data() + static_cast<std::size_t>(y * w + x) * 4);
    }
    if (rOpt.mGridFloor && h - horizon > 3) {
        const int phase = ((rOpt.mThemePhase % 8) + 8) % 8;
        const RndColor line = (phase % 2 == 0) ? RndColor{255, 50, 200, 255} : RndColor{0, 220, 255, 255};
        const double yh = static_cast<double>(horizon);
        const double yb = static_cast<double>(h - 1);
        const int rows = 8;
        for (int k = 0; k < rows; ++k) {
            const double d = (static_cast<double>(k) + static_cast<double>(phase) / 8.0) / rows;
            const double y = yh + (yb - yh) * d * d;
            rnd_line(bg, 0.0, y, static_cast<double>(w - 1), y, line);
        }
        const int rays = 16;
        for (int j = 0; j <= rays; ++j) {
            const double xb = static_cast<double>(w) / 2.0 +
                              (static_cast<double>(j) - rays / 2.0) * static_cast<double>(w) * 2.0 / rays;
            rnd_line(bg, static_cast<double>(w) / 2.0, yh, xb, yb, line);
        }
    }
    for (std::size_t i = 0; i < static_cast<std::size_t>(w * h); ++i) {
        std::uint8_t* p = rFrame.mRgba.data() + i * 4;
        const int a = p[3];
        if (a == 255)
            continue;
        const std::uint8_t* q = bg.mRgba.data() + i * 4;
        for (int c = 0; c < 3; ++c)
            p[c] = static_cast<std::uint8_t>((p[c] * a + q[c] * (255 - a) + 127) / 255);
        p[3] = 255;
    }
}

// Bloom, fringe and scanlines on the final frame, in that order.
void rnd_theme_post(Frame& rFrame, const RenderOptions& rOpt) {
    const std::int64_t w = rFrame.mWidth;
    const std::int64_t h = rFrame.mHeight;
    if (w < 1 || h < 1)
        return;
    std::uint8_t* px = rFrame.mRgba.data();
    if (rOpt.mBloom) {
        const std::int64_t radius = std::max<std::int64_t>(1, std::min(w, h) / 100);
        std::vector<std::int32_t> bright(static_cast<std::size_t>(w * h) * 3, 0);
        for (std::size_t i = 0; i < static_cast<std::size_t>(w * h); ++i) {
            const std::uint8_t* p = px + i * 4;
            const int lum = (p[0] * 54 + p[1] * 183 + p[2] * 19) >> 8;
            if (lum >= 190)
                for (int c = 0; c < 3; ++c)
                    bright[i * 3 + static_cast<std::size_t>(c)] = p[c];
        }
        // A separable box blur: sums over a window of 2 * radius + 1, divided once.
        std::vector<std::int32_t> tmp(bright.size(), 0);
        const std::int64_t window = 2 * radius + 1;
        for (std::int64_t y = 0; y < h; ++y)
            for (int c = 0; c < 3; ++c) {
                std::int64_t sum = 0;
                for (std::int64_t x = -radius; x <= radius; ++x)
                    if (x >= 0 && x < w)
                        sum += bright[static_cast<std::size_t>(y * w + x) * 3 + static_cast<std::size_t>(c)];
                for (std::int64_t x = 0; x < w; ++x) {
                    tmp[static_cast<std::size_t>(y * w + x) * 3 + static_cast<std::size_t>(c)] =
                        static_cast<std::int32_t>(sum / window);
                    const std::int64_t add = x + radius + 1;
                    const std::int64_t drop = x - radius;
                    if (add < w)
                        sum += bright[static_cast<std::size_t>(y * w + add) * 3 + static_cast<std::size_t>(c)];
                    if (drop >= 0)
                        sum -= bright[static_cast<std::size_t>(y * w + drop) * 3 + static_cast<std::size_t>(c)];
                }
            }
        std::vector<std::int32_t> blur(bright.size(), 0);
        for (std::int64_t x = 0; x < w; ++x)
            for (int c = 0; c < 3; ++c) {
                std::int64_t sum = 0;
                for (std::int64_t y = -radius; y <= radius; ++y)
                    if (y >= 0 && y < h)
                        sum += tmp[static_cast<std::size_t>(y * w + x) * 3 + static_cast<std::size_t>(c)];
                for (std::int64_t y = 0; y < h; ++y) {
                    blur[static_cast<std::size_t>(y * w + x) * 3 + static_cast<std::size_t>(c)] =
                        static_cast<std::int32_t>(sum / window);
                    const std::int64_t add = y + radius + 1;
                    const std::int64_t drop = y - radius;
                    if (add < h)
                        sum += tmp[static_cast<std::size_t>(add * w + x) * 3 + static_cast<std::size_t>(c)];
                    if (drop >= 0)
                        sum -= tmp[static_cast<std::size_t>(drop * w + x) * 3 + static_cast<std::size_t>(c)];
                }
            }
        for (std::size_t i = 0; i < static_cast<std::size_t>(w * h); ++i)
            for (int c = 0; c < 3; ++c)
                px[i * 4 + static_cast<std::size_t>(c)] = static_cast<std::uint8_t>(
                    std::min<std::int32_t>(255, px[i * 4 + static_cast<std::size_t>(c)] +
                                                    blur[i * 3 + static_cast<std::size_t>(c)] / 2));
    }
    if (rOpt.mFringe) {
        const std::int64_t d = std::max<std::int64_t>(1, w / 320);
        const std::vector<std::uint8_t> src(px, px + static_cast<std::size_t>(w * h) * 4);
        for (std::int64_t y = 0; y < h; ++y)
            for (std::int64_t x = 0; x < w; ++x) {
                const std::size_t at = static_cast<std::size_t>(y * w + x) * 4;
                const std::size_t left = static_cast<std::size_t>(y * w + std::max<std::int64_t>(0, x - d)) * 4;
                const std::size_t right = static_cast<std::size_t>(y * w + std::min<std::int64_t>(w - 1, x + d)) * 4;
                px[at] = src[left];
                px[at + 2] = src[right + 2];
            }
    }
    if (rOpt.mScanlines)
        for (std::int64_t y = 1; y < h; y += 2)
            for (std::int64_t x = 0; x < w; ++x) {
                std::uint8_t* p = px + static_cast<std::size_t>(y * w + x) * 4;
                for (int c = 0; c < 3; ++c)
                    p[c] = static_cast<std::uint8_t>(p[c] - p[c] / 4);
            }
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

namespace {

struct RndExtra {
    std::int64_t mA;
    std::int64_t mB;
    RndColor mColor;
};

/// One arrow of the vector layer: its barbs depend on the camera, so only the
/// shaft is fixed.
struct RndArrow {
    double mP[3];
    double mTip[3];
    double mDir[3];
    double mLen = 0.0;
};

/// Everything a frame needs that does not depend on the camera: the drawn
/// faces, the skin of a volume, the resolved field and its range, normals,
/// edge lines and extra layers. `render_scene` draws it from any viewpoint.
struct RndPrepared {
    RenderOptions mOpt;  ///< the options it was prepared with (field, edges, colours)
    std::size_t mNumFaces = 0;
    std::vector<std::int64_t> mFStart;
    std::vector<std::int64_t> mFNodes;
    std::vector<std::int64_t> mFIds;
    std::vector<std::int64_t> mSNodes;  ///< the rings through smooth shading's split points
    std::vector<double> mXyz;           ///< points, then the extra layers' vertices
    std::vector<double> mFNormal;
    std::vector<double> mVNormal;
    bool mSmooth = false;
    RndGroup mLines;
    RndGroup mPoints;
    std::vector<double> mFValues;
    std::vector<double> mLValues;
    std::vector<double> mPValues;
    RndColoring mColoring;
    bool mContinuous = false;
    bool mCategorical = false;
    std::string mLabel;
    std::optional<int> mComponent;
    std::vector<std::string> mKeys;
    std::vector<RndExtra> mExtras;
    std::size_t mExtrasBeforeArrows = 0;
    std::vector<RndArrow> mArrows;
    std::vector<detail::RasterLine> mEdgeLines;
    bool mDeferCategoryEdges = false;  ///< orientation categories need the camera
    std::array<double, 3> mCentre = {0.0, 0.0, 0.0};  ///< bounding sphere of what is drawn
    double mRadius = 1.0;
};

/// The lines between faces whose categories differ: every interior edge shared by
/// exactly two faces of different value.
void rnd_category_edge_lines(const RenderOptions& rOpt, bool Categorical,
                             const std::vector<std::int64_t>& rFStart,
                             const std::vector<std::int64_t>& rFNodes,
                             const std::vector<double>& rFValues,
                             std::vector<detail::RasterLine>& rLines) {
    const std::vector<std::int64_t>& f_start = rFStart;
    const std::vector<std::int64_t>& f_nodes = rFNodes;
    const std::vector<double>& f_values = rFValues;
    const std::size_t num_faces = f_start.empty() ? 0 : f_start.size() - 1;
    const bool categorical = Categorical;
    if (!categorical)
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "category edges need a categorical colouring "
                                    "(categorical, color_regions or a flag diagnostic)");
    struct Use {
        std::int64_t mLo, mHi;
        std::size_t mFace;
        bool operator<(const Use& rO) const {
            return mLo != rO.mLo ? mLo < rO.mLo
                                 : (mHi != rO.mHi ? mHi < rO.mHi : mFace < rO.mFace);
        }
    };
    std::vector<Use> uses;
    for (std::size_t f = 0; f < num_faces; ++f) {
        const std::int64_t lo = f_start[f];
        const std::int64_t n = f_start[f + 1] - lo;
        for (std::int64_t k = 0; k < n; ++k) {
            const std::int64_t a = f_nodes[static_cast<std::size_t>(lo + k)];
            const std::int64_t b = f_nodes[static_cast<std::size_t>(lo + (k + 1) % n)];
            if (a != b)
                uses.push_back({std::min(a, b), std::max(a, b), f});
        }
    }
    parallel_sort(uses.begin(), uses.end());
    auto same = [&](double a, double b) {
        return (std::isfinite(a) ? a : -1.0) == (std::isfinite(b) ? b : -1.0);
    };
    for (std::size_t i = 0; i < uses.size();) {
        std::size_t j = i;
        while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
            ++j;
        if (j - i == 2 && !same(f_values[uses[i].mFace], f_values[uses[i + 1].mFace]))
            rLines.push_back({uses[i].mLo, uses[i].mHi, rOpt.mEdgeColor, -1});
        i = j;
    }
}

// ---------------------------------------------------------------------------
// Cut-aways: clip the drawn geometry against up to two planes
// ---------------------------------------------------------------------------

// The kept side of a plane: dist(x) = n . x + d >= 0, n a unit vector.
struct RndPlane {
    double mN[3];
    double mD;
};

std::vector<RndPlane> rnd_planes(const std::vector<RenderCutaway>& rCutaways) {
    std::vector<RndPlane> out;
    for (const RenderCutaway& cut : rCutaways) {
        RndPlane plane;
        const double len = std::sqrt(cut.mNormal[0] * cut.mNormal[0] +
                                     cut.mNormal[1] * cut.mNormal[1] +
                                     cut.mNormal[2] * cut.mNormal[2]);
        for (std::size_t k = 0; k < 3; ++k)
            plane.mN[k] = cut.mNormal[k] / len;
        plane.mD = -(plane.mN[0] * cut.mPoint[0] + plane.mN[1] * cut.mPoint[1] +
                     plane.mN[2] * cut.mPoint[2]);
        out.push_back(plane);
    }
    return out;
}

double rnd_dist(const RndPlane& rPlane, const double* pX) {
    return rPlane.mN[0] * pX[0] + rPlane.mN[1] * pX[1] + rPlane.mN[2] * pX[2] + rPlane.mD;
}

// A polygon corner while it is being clipped: its position, its smooth-shading
// normal, and the vertex it came from (-1 for a point made by the cut).
struct RndClipCorner {
    double mX[3];
    double mN[3];
    std::int64_t mId;
};

// Sutherland-Hodgman against one plane. The new corners lie where an edge
// meets the plane, with the position and the normal interpolated along it.
void rnd_clip_polygon(const RndPlane& rPlane, std::vector<RndClipCorner>& rPoly) {
    std::vector<RndClipCorner> out;
    const std::size_t n = rPoly.size();
    for (std::size_t i = 0; i < n; ++i) {
        const RndClipCorner& a = rPoly[i];
        const RndClipCorner& b = rPoly[(i + 1) % n];
        const double da = rnd_dist(rPlane, a.mX);
        const double db = rnd_dist(rPlane, b.mX);
        if (da >= 0.0)
            out.push_back(a);
        if ((da >= 0.0) != (db >= 0.0)) {
            const double t = da / (da - db);
            RndClipCorner q;
            q.mId = -1;
            for (std::size_t k = 0; k < 3; ++k) {
                q.mX[k] = a.mX[k] + t * (b.mX[k] - a.mX[k]);
                q.mN[k] = a.mN[k] + t * (b.mN[k] - a.mN[k]);
            }
            const double len = std::sqrt(q.mN[0] * q.mN[0] + q.mN[1] * q.mN[1] + q.mN[2] * q.mN[2]);
            if (len > 0.0 && std::isfinite(len))
                for (double& c : q.mN)
                    c /= len;
            out.push_back(q);
        }
    }
    rPoly.swap(out);
}

// The geometry left once the planes have cut it.
struct RndClipped {
    std::size_t mNumFaces = 0;
    std::vector<std::int64_t> mFStart = {0};
    std::vector<std::int64_t> mFNodes;  // the rasterized rings (also used for the diagonals)
    std::vector<std::int64_t> mFIds;
    std::vector<double> mFNormal;
    std::vector<double> mFValues;
    std::vector<double> mVNormal;  // grown with the new corners, NaN where unknown
    struct Line {
        std::int64_t mA;
        std::int64_t mB;
        std::size_t mSource;
    };
    std::vector<Line> mLines;
    std::vector<std::size_t> mPoints;  // the input points that survive
    std::vector<detail::RasterLine> mEdgeLines;
};

// A segment against every plane; false when nothing of it is kept. The ends the
// cut moves become new vertices appended to rXyz.
bool rnd_clip_segment(const std::vector<RndPlane>& rPlanes, std::vector<double>& rXyz,
                      std::int64_t& rA, std::int64_t& rB) {
    double pa[3];
    double pb[3];
    for (std::size_t k = 0; k < 3; ++k) {
        pa[k] = rXyz[3 * static_cast<std::size_t>(rA) + k];
        pb[k] = rXyz[3 * static_cast<std::size_t>(rB) + k];
    }
    bool moved_a = false;
    bool moved_b = false;
    for (const RndPlane& plane : rPlanes) {
        const double da = rnd_dist(plane, pa);
        const double db = rnd_dist(plane, pb);
        if (da < 0.0 && db < 0.0)
            return false;
        if (da >= 0.0 && db >= 0.0)
            continue;
        const double t = da / (da - db);
        double q[3];
        for (std::size_t k = 0; k < 3; ++k)
            q[k] = pa[k] + t * (pb[k] - pa[k]);
        if (da < 0.0) {
            for (std::size_t k = 0; k < 3; ++k)
                pa[k] = q[k];
            moved_a = true;
        } else {
            for (std::size_t k = 0; k < 3; ++k)
                pb[k] = q[k];
            moved_b = true;
        }
    }
    if (moved_a) {
        rA = static_cast<std::int64_t>(rXyz.size() / 3);
        rXyz.insert(rXyz.end(), pa, pa + 3);
    }
    if (moved_b) {
        rB = static_cast<std::int64_t>(rXyz.size() / 3);
        rXyz.insert(rXyz.end(), pb, pb + 3);
    }
    return true;
}

bool rnd_point_kept(const std::vector<RndPlane>& rPlanes, const double* pX) {
    for (const RndPlane& plane : rPlanes)
        if (rnd_dist(plane, pX) < 0.0)
            return false;
    return true;
}

// The options with the theme's default colours in place of the ones left at
// their defaults (a colour you set to the default itself cannot be told apart,
// and takes the theme's).
RenderOptions rnd_themed(const RenderOptions& rIn) {
    RenderOptions out = rIn;
    if (out.mTheme == RenderTheme::Synthwave) {
        const RenderOptions defaults;
        if (out.mFillColor == defaults.mFillColor)
            out.mFillColor = {96, 56, 190, 255};
        if (out.mEdgeColor == defaults.mEdgeColor)
            out.mEdgeColor = {0, 230, 255, 255};
        if (out.mLineColor == defaults.mLineColor)
            out.mLineColor = {255, 60, 200, 255};
    }
    return out;
}

std::shared_ptr<const RndPrepared> rnd_prepare(const Mesh& rMesh, const RenderOptions& rOptIn) {
    const RenderOptions rOpt = rnd_themed(rOptIn);
    rnd_validate(rOpt);
    double azimuth = 0.0;
    double elevation = 0.0;
    detail::render_view_angles(rOpt, azimuth, elevation);  // rejects an unknown view

    // The array that is coloured: an input array, or one derived from it.
    Mesh work;
    const Mesh* p_src = &rMesh;
    std::string color_by;
    std::string label;
    std::optional<int> component;
    rnd_prepare_source(rOpt, rMesh, work, p_src, color_by, label, component);
    const Mesh& src = *p_src;
    if (rOpt.mCategorical && color_by.empty())
        throw std::invalid_argument(std::string(kRndPrefix) + "categorical needs color_by");
    if ((rOpt.mIsolines > 0 || !rOpt.mIsoLevels.empty()) && color_by.empty())
        throw std::invalid_argument(std::string(kRndPrefix) +
                                    "isolines need color_by (a point array)");

    Mesh skin;
    bool has_skin = false;
    RndGeometry g = rnd_gather(src, skin, has_skin);
    const std::size_t num_source_points = g.mNumSource;
    const Mesh* p_skin = has_skin ? &skin : nullptr;

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

    // A warp moves the points before anything is measured from them; the
    // undeformed positions are kept for the outline.
    std::vector<double> xyz_undeformed;
    if (!rOpt.mWarp.empty()) {
        xyz_undeformed = g.mXyz;
        const std::vector<double> disp =
            rnd_vertex_array(src, p_skin, num_source_points, rOpt.mWarp, 3, "warp");
        for (std::size_t i = 0; i < g.mXyz.size(); ++i)
            g.mXyz[i] += rOpt.mWarpScale * disp[i];
    }

    // Face normals (Newell's, from the shared crease kernel).
    std::vector<double> f_normal(3 * num_faces, 0.0);
    parallel_for(num_faces, [&](std::size_t f) {
        detail::ring_unit_normal(g.mXyz.data(), f_nodes.data() + f_start[f],
                                 static_cast<std::size_t>(f_start[f + 1] - f_start[f]),
                                 &f_normal[3 * f]);
    });

    // --- Values -----------------------------------------------------------
    RndColoring coloring;
    coloring.mNan = rOpt.mNanColor;
    coloring.mScale = rOpt.mScale;
    coloring.mThreshold = rOpt.mScaleThreshold;
    std::vector<double> f_values(num_faces, 0.0);
    std::vector<double> l_values(g.mLines.mIds.size(), 0.0);
    std::vector<double> p_values(g.mPoints.mIds.size(), 0.0);
    std::vector<std::string> keys;
    bool continuous = false;
    bool categorical = false;
    const bool diag_flag = rOpt.mDiagnostic == RenderDiagnostic::Inverted ||
                           rOpt.mDiagnostic == RenderDiagnostic::Degenerate;
    detail::ColorSpec spec;
    spec.mColorBy = color_by;
    spec.mComponent = component;
    spec.mCmap = rOpt.mCmap;

    // The first Cell region that holds each input cell.
    std::vector<std::int32_t> cell_region;
    std::vector<std::string> region_names;
    std::vector<std::size_t> region_counts;

    if (rOpt.mColorRegions) {
        std::size_t total_cells = 0;
        for (const auto cb : src.CellRange())
            total_cells += cb.NumCells();
        cell_region.assign(total_cells, -1);
        for (std::size_t r = 0; r < src.NumRegions(); ++r) {
            const meshioplusplus::Region& region = src.Region(r);
            if (region.mKind != RegionKind::Cell)
                continue;
            const std::int32_t index = static_cast<std::int32_t>(region_names.size());
            region_names.push_back(region.mName);
            region_counts.push_back(region.NumEntries());
            for (std::size_t e = 0; e < region.mEntries.Size(); ++e) {
                const std::int64_t c = detail::read_int(region.mEntries, e);
                if (c >= 0 && static_cast<std::size_t>(c) < total_cells &&
                    cell_region[static_cast<std::size_t>(c)] < 0)
                    cell_region[static_cast<std::size_t>(c)] = index;
            }
        }
        if (region_names.empty())
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "color_regions needs at least one named cell region");
        auto region_of = [&](std::int64_t id) {
            return id >= 0 && static_cast<std::size_t>(id) < cell_region.size() &&
                           cell_region[static_cast<std::size_t>(id)] >= 0
                       ? static_cast<double>(cell_region[static_cast<std::size_t>(id)])
                       : std::nan("");
        };
        for (std::size_t f = 0; f < num_faces; ++f)
            f_values[f] = region_of(f_ids[f]);
        for (std::size_t i = 0; i < l_values.size(); ++i)
            l_values[i] = region_of(g.mLines.mIds[i]);
        for (std::size_t i = 0; i < p_values.size(); ++i)
            p_values[i] = region_of(g.mPoints.mIds[i]);
        categorical = true;
        for (std::size_t r = 0; r < region_names.size(); ++r) {
            const auto& pal = detail::render_category_palette();
            const RndColor c = pal[r % pal.size()];
            coloring.mPalette.push_back(c);
            keys.push_back("region " + region_names[r] + ": " + rnd_hex(c) + " (" +
                           std::to_string(region_counts[r]) + " cells)");
        }
    } else if (!color_by.empty()) {
        coloring.mpTable = detail::colormap_table(rOpt.mCmap);
        const bool is_point = src.HasPointData(color_by);
        if (!is_point && !src.HasCellData(color_by))
            (void)detail::resolve_face_colors(spec, src, src, {});  // throws the shared message
        f_values = rnd_group_values(spec, src, g.mSkinFaces, is_point);
        const std::vector<double> own = rnd_group_values(spec, src, g.mFaces, is_point);
        f_values.insert(f_values.end(), own.begin(), own.end());
        l_values = rnd_group_values(spec, src, g.mLines, is_point);
        p_values = rnd_group_values(spec, src, g.mPoints, is_point);
        if (diag_flag) {
            // 0 or 1; a cell the metric does not apply to is not flagged.
            for (std::vector<double>* p_list : {&f_values, &l_values, &p_values})
                for (double& v : *p_list)
                    v = std::isfinite(v) && v > 0.5 ? 1.0 : 0.0;
            coloring.mPalette = {rOpt.mFillColor, {220, 40, 40, 255}};
            coloring.mFacesOnly = true;
            std::size_t flagged = 0;
            for (double v : f_values)
                flagged += v > 0.5 ? 1 : 0;
            keys.push_back(std::string(rOpt.mDiagnostic == RenderDiagnostic::Inverted
                                           ? "inverted"
                                           : "degenerate") +
                           ": " + std::to_string(flagged) + " drawn faces in red");
            categorical = true;
        } else if (rOpt.mCategorical) {
            std::vector<double> distinct;
            for (const std::vector<double>* p_list : {&f_values, &l_values, &p_values})
                for (double v : *p_list)
                    if (std::isfinite(v))
                        distinct.push_back(v);
            std::sort(distinct.begin(), distinct.end());
            distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
            std::vector<std::size_t> counts(distinct.size(), 0);
            auto index_of = [&](double v) {
                return static_cast<std::size_t>(
                    std::lower_bound(distinct.begin(), distinct.end(), v) - distinct.begin());
            };
            for (double v : f_values)
                if (std::isfinite(v))
                    ++counts[index_of(v)];
            for (std::vector<double>* p_list : {&f_values, &l_values, &p_values})
                for (double& v : *p_list)
                    if (std::isfinite(v))
                        v = static_cast<double>(index_of(v));
            const auto& pal = detail::render_category_palette();
            for (std::size_t i = 0; i < distinct.size(); ++i) {
                const RndColor c = pal[i % pal.size()];
                coloring.mPalette.push_back(c);
                keys.push_back(label + " = " + rnd_num(distinct[i]) + ": " + rnd_hex(c) + " (" +
                               std::to_string(counts[i]) + " drawn faces)");
            }
            if (coloring.mPalette.empty())
                coloring.mPalette.push_back(rOpt.mFillColor);
            categorical = true;
        } else {
            continuous = true;
        }
    } else if (rOpt.mDiagnostic == RenderDiagnostic::EdgeLength) {
        coloring.mpTable = detail::colormap_table(rOpt.mCmap);
        label = "edge length";
        for (std::size_t f = 0; f < num_faces; ++f) {
            const std::int64_t lo = f_start[f];
            const std::int64_t n = f_start[f + 1] - lo;
            double sum = 0.0;
            for (std::int64_t k = 0; k < n; ++k) {
                const std::size_t a =
                    static_cast<std::size_t>(f_nodes[static_cast<std::size_t>(lo + k)]);
                const std::size_t b =
                    static_cast<std::size_t>(f_nodes[static_cast<std::size_t>(lo + (k + 1) % n)]);
                double d2 = 0.0;
                for (std::size_t c = 0; c < 3; ++c) {
                    const double d = g.mXyz[3 * a + c] - g.mXyz[3 * b + c];
                    d2 += d * d;
                }
                sum += std::sqrt(d2);
            }
            f_values[f] = n > 0 ? sum / static_cast<double>(n) : std::nan("");
        }
        std::fill(l_values.begin(), l_values.end(), std::nan(""));
        std::fill(p_values.begin(), p_values.end(), std::nan(""));
        coloring.mFacesOnly = true;
        continuous = true;
    }

    if (continuous) {
        // The range: the finite values of what is drawn (positive ones on a log
        // scale), narrowed to percentiles, made symmetric, or given outright.
        std::vector<double> finite;
        for (const std::vector<double>* p_list : {&f_values, &l_values, &p_values})
            for (double v : *p_list)
                if (std::isfinite(v) && (rOpt.mScale != RenderScale::Log || v > 0.0))
                    finite.push_back(v);
        double lo = 0.0;
        double hi = 0.0;
        if (!finite.empty()) {
            if (rOpt.mClipLow.has_value() || rOpt.mClipHigh.has_value()) {
                std::vector<double> sorted = finite;
                lo = detail::render_percentile(sorted, rOpt.mClipLow.value_or(0.0));
                hi = detail::render_percentile(sorted, rOpt.mClipHigh.value_or(100.0));
            } else {
                const auto mm = std::minmax_element(finite.begin(), finite.end());
                lo = *mm.first;
                hi = *mm.second;
            }
        }
        if (rOpt.mSymmetric) {
            const double m = std::max(std::fabs(lo), std::fabs(hi));
            lo = -m;
            hi = m;
        }
        coloring.mVMin = rOpt.mVMin.value_or(lo);
        coloring.mVMax = rOpt.mVMax.value_or(hi);
        if (coloring.mVMin > coloring.mVMax)
            throw std::invalid_argument(std::string(kRndPrefix) + "vmin must not exceed vmax");
        if (rOpt.mScale == RenderScale::Log && !(coloring.mVMin > 0.0) && !finite.empty())
            throw std::invalid_argument(std::string(kRndPrefix) +
                                        "a log scale needs a positive range, not " +
                                        rnd_num(coloring.mVMin) + " .. " + rnd_num(coloring.mVMax));
    }
    coloring.mActive = continuous || categorical;

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

    // --- Extra line layers: contours, arrows, the undeformed outline ------
    struct Extra {
        std::int64_t mA;
        std::int64_t mB;
        RndColor mColor;
    };
    std::vector<Extra> extras;
    auto add_vertex = [&](const double* pPoint) {
        const std::int64_t id = static_cast<std::int64_t>(xyz.size() / 3);
        xyz.insert(xyz.end(), pPoint, pPoint + 3);
        return id;
    };
    auto add_segment = [&](const double* pA, const double* pB, const RndColor& rColor) {
        const std::int64_t a = add_vertex(pA);
        const std::int64_t b = add_vertex(pB);
        extras.push_back({a, b, rColor});
    };

    if (rOpt.mIsolines > 0 || !rOpt.mIsoLevels.empty()) {
        if (!src.HasPointData(color_by))
            throw std::invalid_argument(
                std::string(kRndPrefix) + "isolines need a point array, but '" + color_by +
                "' is not point data (available: " + rnd_available(src) + ")");
        const std::vector<double> vv = rnd_vertex_scalar(spec, src, p_skin, num_source_points);
        std::vector<double> levels = rOpt.mIsoLevels;
        if (levels.empty()) {
            if (!continuous)
                throw std::invalid_argument(std::string(kRndPrefix) +
                                            "isolines need a continuous range to divide");
            for (std::int32_t k = 1; k <= rOpt.mIsolines; ++k)
                levels.push_back(coloring.mVMin + (coloring.mVMax - coloring.mVMin) *
                                                      static_cast<double>(k) /
                                                      static_cast<double>(rOpt.mIsolines + 1));
        }
        std::vector<std::array<int, 3>> tris;
        for (std::size_t f = 0; f < num_faces; ++f) {
            const std::int64_t* ring = f_nodes.data() + f_start[f];
            rnd_face_triangles(xyz, ring, f_start[f + 1] - f_start[f], tris);
            for (const std::array<int, 3>& t : tris) {
                double corner[9];
                double value[3];
                for (int k = 0; k < 3; ++k) {
                    const std::size_t v =
                        static_cast<std::size_t>(ring[t[static_cast<std::size_t>(k)]]);
                    for (std::size_t c = 0; c < 3; ++c)
                        corner[3 * k + static_cast<int>(c)] = xyz[3 * v + c];
                    value[k] = vv[v];
                }
                for (double level : levels) {
                    double seg[6];
                    if (detail::render_triangle_contour(corner, value, level, seg))
                        add_segment(seg, seg + 3, rOpt.mIsoColor);
                }
            }
        }
        keys.push_back("isolines: " + std::to_string(levels.size()) + " levels (" +
                       rnd_hex(rOpt.mIsoColor) + ")");
    }

    if (!rOpt.mStreamlines.empty()) {
        const std::vector<double> vec =
            rnd_vertex_array(src, p_skin, num_source_points, rOpt.mStreamlines, 3, "streamline");
        detail::StreamlineOptions stream;
        stream.mSeeds = static_cast<std::size_t>(rOpt.mStreamSeeds);
        stream.mLength = rOpt.mStreamLength;
        detail::Streamlines lines;
        try {
            lines = detail::trace_streamlines(src, xyz.data(), vec.data(), stream);
        } catch (const std::invalid_argument& rErr) {
            throw std::invalid_argument(std::string(kRndPrefix) + rErr.what());
        }
        for (std::size_t l = 0; l < lines.NumLines(); ++l)
            for (std::size_t k = lines.mStart[l]; k + 1 < lines.mStart[l + 1]; ++k)
                add_segment(&lines.mXyz[3 * k], &lines.mXyz[3 * (k + 1)], rOpt.mStreamColor);
        keys.push_back("streamlines: " + rOpt.mStreamlines + ", " +
                       std::to_string(lines.NumLines()) + " lines (" + rnd_hex(rOpt.mStreamColor) +
                       ")");
    }

    const std::size_t extras_before_arrows = extras.size();
    std::vector<RndArrow> arrows;
    if (!rOpt.mVectors.empty() && num_faces > 0) {
        const std::vector<double> vec =
            rnd_vertex_array(src, p_skin, num_source_points, rOpt.mVectors, 3, "vector");
        std::vector<std::int64_t> ids(f_nodes);
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        const std::size_t want =
            std::min<std::size_t>(ids.size(), static_cast<std::size_t>(rOpt.mVectorCount));
        std::vector<std::int64_t> picked;
        for (std::size_t k = 0; k < want; ++k)
            picked.push_back(ids[static_cast<std::size_t>((static_cast<double>(k) + 0.5) *
                                                          static_cast<double>(ids.size()) /
                                                          static_cast<double>(want))]);
        double lo3[3] = {0, 0, 0};
        double hi3[3] = {0, 0, 0};
        bool first = true;
        for (std::int64_t v : ids)
            for (std::size_t c = 0; c < 3; ++c) {
                const double x = xyz[3 * static_cast<std::size_t>(v) + c];
                lo3[c] = first ? x : std::min(lo3[c], x);
                hi3[c] = first ? x : std::max(hi3[c], x);
                first = false;
            }
        first = true;
        const double diag = std::sqrt((hi3[0] - lo3[0]) * (hi3[0] - lo3[0]) +
                                      (hi3[1] - lo3[1]) * (hi3[1] - lo3[1]) +
                                      (hi3[2] - lo3[2]) * (hi3[2] - lo3[2]));
        double max_mag = 0.0;
        for (std::int64_t v : picked) {
            const double* a = &vec[3 * static_cast<std::size_t>(v)];
            const double m = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
            if (std::isfinite(m))
                max_mag = std::max(max_mag, m);
        }
        for (std::int64_t v : picked) {
            const double* a = &vec[3 * static_cast<std::size_t>(v)];
            const double mag = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
            if (!(mag > 0.0) || !std::isfinite(mag) || !(max_mag > 0.0))
                continue;
            const double len =
                rOpt.mVectorLength > 0.0 ? rOpt.mVectorLength : 0.06 * diag * mag / max_mag;
            const double dir[3] = {a[0] / mag, a[1] / mag, a[2] / mag};
            const double* p = &xyz[3 * static_cast<std::size_t>(v)];
            const double tip[3] = {p[0] + dir[0] * len, p[1] + dir[1] * len, p[2] + dir[2] * len};
            RndArrow arrow;
            for (std::size_t c = 0; c < 3; ++c) {
                arrow.mP[c] = p[c];
                arrow.mTip[c] = tip[c];
                arrow.mDir[c] = dir[c];
            }
            arrow.mLen = len;
            arrows.push_back(arrow);
        }
        keys.push_back("vectors: " + rOpt.mVectors + ", " + std::to_string(picked.size()) +
                       " arrows, longest " + rnd_num(max_mag));
    }

    if (!rOpt.mWarp.empty() && rOpt.mWarpOutline && num_faces > 0) {
        std::vector<double> n0(3 * num_faces, 0.0);
        for (std::size_t f = 0; f < num_faces; ++f)
            detail::ring_unit_normal(xyz_undeformed.data(), f_nodes.data() + f_start[f],
                                     static_cast<std::size_t>(f_start[f + 1] - f_start[f]),
                                     &n0[3 * f]);
        for (const detail::CreaseEdge& e :
             detail::crease_edges(f_start, f_nodes, n0, rOpt.mFeatureAngle))
            if (e.IsCrease() || e.IsBoundary())
                add_segment(&xyz_undeformed[3 * static_cast<std::size_t>(e.mLo)],
                            &xyz_undeformed[3 * static_cast<std::size_t>(e.mHi)],
                            rOpt.mOutlineColor);
        keys.push_back("warp: " + rOpt.mWarp + " x " + rnd_num(rOpt.mWarpScale) +
                       ", undeformed outline in " + rnd_hex(rOpt.mOutlineColor));
    } else if (!rOpt.mWarp.empty()) {
        keys.push_back("warp: " + rOpt.mWarp + " x " + rnd_num(rOpt.mWarpScale));
    }

    // Edge overlays that depend on the geometry alone.
    std::vector<detail::RasterLine> edge_lines;
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
            edge_lines.push_back({e.first, e.second, rOpt.mEdgeColor, -1});
    } else if (rOpt.mEdges == RenderEdges::Feature && num_faces > 0) {
        const std::vector<detail::CreaseEdge> edges =
            detail::crease_edges(f_start, f_nodes, f_normal, rOpt.mFeatureAngle);
        for (const detail::CreaseEdge& e : edges)
            if (e.IsCrease() || e.IsBoundary())
                edge_lines.push_back({e.mLo, e.mHi, rOpt.mEdgeColor, -1});
    }
    if (rOpt.mDiagnostic == RenderDiagnostic::FreeEdges && num_faces > 0) {
        const RndColor open_color = {235, 140, 0, 255};
        const RndColor non_manifold_color = {220, 30, 30, 255};
        const RndColor inconsistent_color = {200, 0, 200, 255};
        std::size_t open_edges = 0, non_manifold = 0, inconsistent = 0;
        for (const detail::CreaseEdge& e :
             detail::crease_edges(f_start, f_nodes, f_normal, 180.0)) {
            if (e.IsNonManifold()) {
                edge_lines.push_back({e.mLo, e.mHi, non_manifold_color, -1});
                ++non_manifold;
            } else if (e.IsBoundary()) {
                edge_lines.push_back({e.mLo, e.mHi, open_color, -1});
                ++open_edges;
            } else if (e.mInconsistent) {
                edge_lines.push_back({e.mLo, e.mHi, inconsistent_color, -1});
                ++inconsistent;
            }
        }
        keys.push_back(
            "free edges: " + std::to_string(open_edges) + " open " + rnd_hex(open_color) + ", " +
            std::to_string(non_manifold) + " non-manifold " + rnd_hex(non_manifold_color) + ", " +
            std::to_string(inconsistent) + " inconsistent " + rnd_hex(inconsistent_color));
    }
    bool defer_category_edges = false;
    if (rOpt.mCategoryEdges) {
        if (rOpt.mDiagnostic == RenderDiagnostic::Orientation)
            defer_category_edges = true;  // its categories depend on the camera
        else
            rnd_category_edge_lines(rOpt, categorical, f_start, f_nodes, f_values, edge_lines);
    }


    // The bounding sphere of what is drawn, for a fixed fit.
    const std::size_t nv = xyz.size() / 3;
    std::vector<std::uint8_t> used(nv, 0);
    for (std::int64_t v : s_nodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : g.mLines.mNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : g.mPoints.mNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (const Extra& e : extras) {
        used[static_cast<std::size_t>(e.mA)] = 1;
        used[static_cast<std::size_t>(e.mB)] = 1;
    }
    std::array<double, 3> lo_b = {0.0, 0.0, 0.0};
    std::array<double, 3> hi_b = {0.0, 0.0, 0.0};
    bool seen = false;
    for (std::size_t i = 0; i < nv; ++i) {
        if (!used[i])
            continue;
        for (std::size_t k = 0; k < 3; ++k) {
            const double c = xyz[3 * i + k];
            lo_b[k] = seen ? std::min(lo_b[k], c) : c;
            hi_b[k] = seen ? std::max(hi_b[k], c) : c;
        }
        seen = true;
    }
    auto prep = std::make_shared<RndPrepared>();
    for (std::size_t k = 0; k < 3; ++k)
        prep->mCentre[k] = 0.5 * (lo_b[k] + hi_b[k]);
    double radius_b = 0.0;
    for (std::size_t i = 0; i < nv; ++i) {
        if (!used[i])
            continue;
        double d2 = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            const double d = xyz[3 * i + k] - prep->mCentre[k];
            d2 += d * d;
        }
        radius_b = std::max(radius_b, std::sqrt(d2));
    }
    prep->mRadius = radius_b > 0.0 ? radius_b : 1.0;

    prep->mOpt = rOpt;
    prep->mNumFaces = num_faces;
    prep->mFStart = std::move(f_start);
    prep->mFNodes = std::move(f_nodes);
    prep->mFIds = std::move(f_ids);
    prep->mSNodes = std::move(s_nodes);
    prep->mXyz = std::move(xyz);
    prep->mFNormal = std::move(f_normal);
    prep->mVNormal = std::move(v_normal);
    prep->mSmooth = smooth;
    prep->mLines = g.mLines;
    prep->mLines.mpDraw = nullptr;
    prep->mPoints = g.mPoints;
    prep->mPoints.mpDraw = nullptr;
    prep->mFValues = std::move(f_values);
    prep->mLValues = std::move(l_values);
    prep->mPValues = std::move(p_values);
    prep->mColoring = coloring;
    prep->mContinuous = continuous;
    prep->mCategorical = categorical;
    prep->mLabel = label;
    prep->mComponent = component;
    prep->mKeys = std::move(keys);
    for (const Extra& e : extras)
        prep->mExtras.push_back({e.mA, e.mB, e.mColor});
    prep->mExtrasBeforeArrows = extras_before_arrows;
    prep->mArrows = std::move(arrows);
    prep->mEdgeLines = std::move(edge_lines);
    prep->mDeferCategoryEdges = defer_category_edges;
    return prep;
}

Frame rnd_draw(const RndPrepared& rPrep, const RenderOptions& rCamera, bool FixedFit) {
    // The field, edge and colour options are the prepared ones; the camera and
    // the frame are the caller's.
    RenderOptions rOpt = rPrep.mOpt;
    rOpt.mWidth = rCamera.mWidth;
    rOpt.mHeight = rCamera.mHeight;
    rOpt.mPixelAspect = rCamera.mPixelAspect;
    rOpt.mSupersample = rCamera.mSupersample;
    rOpt.mAzimuth = rCamera.mAzimuth;
    rOpt.mElevation = rCamera.mElevation;
    rOpt.mRoll = rCamera.mRoll;
    rOpt.mView = rCamera.mView;
    rOpt.mProjection = rCamera.mProjection;
    rOpt.mFovDeg = rCamera.mFovDeg;
    rOpt.mZoom = rCamera.mZoom;
    rOpt.mPanX = rCamera.mPanX;
    rOpt.mPanY = rCamera.mPanY;
    rOpt.mShading = rCamera.mShading;
    rOpt.mTwoSided = rCamera.mTwoSided;
    rOpt.mAmbient = rCamera.mAmbient;
    rOpt.mLightDir = rCamera.mLightDir;
    rOpt.mBackground = rCamera.mBackground;
    rOpt.mPointRadius = rCamera.mPointRadius;
    rOpt.mAxes = rCamera.mAxes;
    rOpt.mScaleBar = rCamera.mScaleBar;
    rOpt.mColorbar = rCamera.mColorbar;
    rOpt.mCutaways = rCamera.mCutaways;
    rOpt.mCutawayTint = rCamera.mCutawayTint;
    rOpt.mTheme = rCamera.mTheme;
    rOpt.mScanlines = rCamera.mScanlines;
    rOpt.mBloom = rCamera.mBloom;
    rOpt.mFringe = rCamera.mFringe;
    rOpt.mGridFloor = rCamera.mGridFloor;
    rOpt.mThemePhase = rCamera.mThemePhase;
    rnd_validate(rOpt);
    double azimuth = 0.0;
    double elevation = 0.0;
    detail::render_view_angles(rOpt, azimuth, elevation);
    const detail::CameraBasis cam = detail::camera_basis(azimuth, elevation, rOpt.mRoll);

    const bool smooth = rPrep.mSmooth && rOpt.mShading == RenderShading::Smooth;
    const std::vector<double>& l_values = rPrep.mLValues;
    const std::vector<double>& p_values = rPrep.mPValues;
    const bool continuous = rPrep.mContinuous;
    bool categorical = rPrep.mCategorical;
    const std::string& label = rPrep.mLabel;
    const std::optional<int>& component = rPrep.mComponent;
    RndColoring coloring = rPrep.mColoring;
    std::vector<std::string> keys = rPrep.mKeys;
    const bool orientation = rOpt.mDiagnostic == RenderDiagnostic::Orientation;

    // The extra layers: the prepared ones with the arrows between them, whose
    // barbs lie in the plane of the arrow and the line of sight.
    std::vector<double> xyz = rPrep.mXyz;
    struct Extra {
        std::int64_t mA;
        std::int64_t mB;
        RndColor mColor;
    };
    std::vector<Extra> extras;
    auto add_vertex = [&](const double* pPoint) {
        const std::int64_t id = static_cast<std::int64_t>(xyz.size() / 3);
        xyz.insert(xyz.end(), pPoint, pPoint + 3);
        return id;
    };
    auto add_segment = [&](const double* pA, const double* pB, const RndColor& rColor) {
        const std::int64_t a = add_vertex(pA);
        const std::int64_t b = add_vertex(pB);
        extras.push_back({a, b, rColor});
    };
    for (std::size_t i = 0; i < rPrep.mExtrasBeforeArrows; ++i)
        extras.push_back({rPrep.mExtras[i].mA, rPrep.mExtras[i].mB, rPrep.mExtras[i].mColor});
    const double cos_barb = 0.9063077870366499;  // 25 degrees
    const double sin_barb = 0.42261826174069944;
    for (const RndArrow& arrow : rPrep.mArrows) {
        const double* dir = arrow.mDir;
        const double len = arrow.mLen;
        double perp[3] = {dir[1] * cam.mW[2] - dir[2] * cam.mW[1],
                          dir[2] * cam.mW[0] - dir[0] * cam.mW[2],
                          dir[0] * cam.mW[1] - dir[1] * cam.mW[0]};
        double pl = std::sqrt(perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]);
        if (pl < 1e-9) {
            perp[0] = cam.mU[0];
            perp[1] = cam.mU[1];
            perp[2] = cam.mU[2];
            pl = 1.0;
        }
        for (double& c : perp)
            c /= pl;
        add_segment(arrow.mP, arrow.mTip, rOpt.mVectorColor);
        for (double sign : {1.0, -1.0}) {
            double barb[3];
            for (std::size_t c = 0; c < 3; ++c)
                barb[c] = arrow.mTip[c] - 0.3 * len * (cos_barb * dir[c] - sign * sin_barb * perp[c]);
            add_segment(arrow.mTip, barb, rOpt.mVectorColor);
        }
    }
    for (std::size_t i = rPrep.mExtrasBeforeArrows; i < rPrep.mExtras.size(); ++i)
        extras.push_back({rPrep.mExtras[i].mA, rPrep.mExtras[i].mB, rPrep.mExtras[i].mColor});

    // Cut-aways: clip the faces, lines and points against the planes. A face the
    // cut leaves whole keeps its vertices; one it crosses is replaced by the
    // polygon on the kept side, with the interpolated corners appended to xyz.
    const std::vector<RndPlane> planes = rnd_planes(rOpt.mCutaways);
    const bool clip_on = !planes.empty();
    RndClipped clipped;
    std::vector<Extra> extras_unclipped;  // the framing is the whole model's, cut or not
    if (clip_on) {
        const double nan = std::nan("");
        if (smooth)
            clipped.mVNormal = rPrep.mVNormal;
        auto new_vertex = [&](const RndClipCorner& rCorner) {
            const std::int64_t id = static_cast<std::int64_t>(xyz.size() / 3);
            xyz.insert(xyz.end(), rCorner.mX, rCorner.mX + 3);
            if (smooth) {
                clipped.mVNormal.resize(3 * static_cast<std::size_t>(id), nan);
                clipped.mVNormal.insert(clipped.mVNormal.end(), rCorner.mN, rCorner.mN + 3);
            }
            return id;
        };
        for (std::size_t f = 0; f < rPrep.mNumFaces; ++f) {
            const std::int64_t lo = rPrep.mFStart[f];
            const std::int64_t n = rPrep.mFStart[f + 1] - lo;
            if (n < 3)
                continue;
            const std::int64_t* ring = rPrep.mSNodes.data() + lo;
            bool inside = true;
            for (std::int64_t k = 0; k < n && inside; ++k)
                for (const RndPlane& plane : planes)
                    if (rnd_dist(plane, &xyz[3 * static_cast<std::size_t>(ring[k])]) < 0.0) {
                        inside = false;
                        break;
                    }
            if (inside) {
                clipped.mFNodes.insert(clipped.mFNodes.end(), ring, ring + n);
            } else {
                std::vector<RndClipCorner> poly;
                for (std::int64_t k = 0; k < n; ++k) {
                    RndClipCorner corner;
                    const std::size_t v = static_cast<std::size_t>(ring[k]);
                    for (std::size_t c = 0; c < 3; ++c) {
                        corner.mX[c] = xyz[3 * v + c];
                        corner.mN[c] = smooth && 3 * v + c < rPrep.mVNormal.size()
                                           ? rPrep.mVNormal[3 * v + c]
                                           : nan;
                    }
                    corner.mId = ring[k];
                    poly.push_back(corner);
                }
                for (const RndPlane& plane : planes)
                    rnd_clip_polygon(plane, poly);
                if (poly.size() < 3)
                    continue;
                for (const RndClipCorner& corner : poly)
                    clipped.mFNodes.push_back(corner.mId >= 0 ? corner.mId : new_vertex(corner));
            }
            clipped.mFStart.push_back(static_cast<std::int64_t>(clipped.mFNodes.size()));
            clipped.mFIds.push_back(rPrep.mFIds[f]);
            clipped.mFNormal.insert(clipped.mFNormal.end(), rPrep.mFNormal.begin() + 3 * f,
                                    rPrep.mFNormal.begin() + 3 * f + 3);
            clipped.mFValues.push_back(rPrep.mFValues[f]);
        }
        clipped.mNumFaces = clipped.mFIds.size();
        for (std::size_t i = 0; i < rPrep.mLines.mIds.size(); ++i) {
            std::int64_t a = rPrep.mLines.mNodes[static_cast<std::size_t>(rPrep.mLines.mStart[i])];
            std::int64_t b =
                rPrep.mLines.mNodes[static_cast<std::size_t>(rPrep.mLines.mStart[i]) + 1];
            if (rnd_clip_segment(planes, xyz, a, b))
                clipped.mLines.push_back({a, b, i});
        }
        for (const detail::RasterLine& edge : rPrep.mEdgeLines) {
            detail::RasterLine line = edge;
            if (rnd_clip_segment(planes, xyz, line.mA, line.mB))
                clipped.mEdgeLines.push_back(line);
        }
        extras_unclipped = extras;
        std::vector<Extra> kept_extras;
        for (const Extra& e : extras) {
            Extra moved = e;
            if (rnd_clip_segment(planes, xyz, moved.mA, moved.mB))
                kept_extras.push_back(moved);
        }
        extras.swap(kept_extras);
        for (std::size_t i = 0; i < rPrep.mPoints.mIds.size(); ++i) {
            const std::size_t v =
                static_cast<std::size_t>(rPrep.mPoints.mNodes[static_cast<std::size_t>(rPrep.mPoints.mStart[i])]);
            if (rnd_point_kept(planes, &xyz[3 * v]))
                clipped.mPoints.push_back(i);
        }
        keys.push_back("cutaway: " + std::to_string(planes.size()) + " plane" +
                       (planes.size() == 1 ? "" : "s") + ", inside in " +
                       rnd_hex(rOpt.mCutawayTint));
    }
    const std::size_t num_faces = clip_on ? clipped.mNumFaces : rPrep.mNumFaces;
    const std::vector<std::int64_t>& f_start = clip_on ? clipped.mFStart : rPrep.mFStart;
    const std::vector<std::int64_t>& f_nodes = clip_on ? clipped.mFNodes : rPrep.mFNodes;
    const std::vector<std::int64_t>& f_ids = clip_on ? clipped.mFIds : rPrep.mFIds;
    const std::vector<std::int64_t>& s_nodes = clip_on ? clipped.mFNodes : rPrep.mSNodes;
    const std::vector<double>& f_normal = clip_on ? clipped.mFNormal : rPrep.mFNormal;
    const std::vector<double>& v_normal = clip_on ? clipped.mVNormal : rPrep.mVNormal;
    const std::vector<double>& face_values = clip_on ? clipped.mFValues : rPrep.mFValues;
    std::vector<double> f_orient;
    if (orientation)
        f_orient = face_values;
    const std::vector<double>& f_values = orientation ? f_orient : face_values;

    // Project every vertex a primitive uses. The frame is fitted to the geometry
    // as it was before any cut (the corners a cut adds lie inside it), so a
    // cut-away keeps the size and place of the whole model instead of enlarging
    // what is left.
    const std::size_t nv = xyz.size() / 3;
    std::vector<std::uint8_t> used(nv, 0);
    for (std::int64_t v : rPrep.mSNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : rPrep.mLines.mNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (std::int64_t v : rPrep.mPoints.mNodes)
        used[static_cast<std::size_t>(v)] = 1;
    for (const Extra& e : clip_on ? extras_unclipped : extras) {
        used[static_cast<std::size_t>(e.mA)] = 1;
        used[static_cast<std::size_t>(e.mB)] = 1;
    }

    const bool perspective = rOpt.mProjection == RenderProjection::Perspective;
    double fit_distance = 0.0;
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
        fit_distance = distance;
        for (int k = 0; k < 3; ++k)
            eye[k] = centre[k] + cam.mW[static_cast<std::size_t>(k)] * distance;
    }

    // Front and back faces, which need the camera.
    if (orientation) {
        const RndColor front = {70, 130, 230, 255};
        const RndColor back = {230, 140, 40, 255};
        coloring.mPalette = {front, back};
        coloring.mFacesOnly = true;
        coloring.mActive = true;
        categorical = true;
        std::size_t back_faces = 0;
        for (std::size_t f = 0; f < num_faces; ++f) {
            const double* n = &f_normal[3 * f];
            double toward[3] = {cam.mW[0], cam.mW[1], cam.mW[2]};
            if (perspective) {
                double c[3] = {0, 0, 0};
                const std::int64_t cnt = f_start[f + 1] - f_start[f];
                for (std::int64_t k = f_start[f]; k < f_start[f + 1]; ++k)
                    for (std::size_t a = 0; a < 3; ++a)
                        c[a] +=
                            xyz[3 * static_cast<std::size_t>(f_nodes[static_cast<std::size_t>(k)]) +
                                a];
                for (std::size_t a = 0; a < 3; ++a)
                    toward[a] = eye[a] - (cnt > 0 ? c[a] / static_cast<double>(cnt) : 0.0);
            }
            const double d = n[0] * toward[0] + n[1] * toward[1] + n[2] * toward[2];
            f_orient[f] = d < 0.0 ? 1.0 : 0.0;
            back_faces += d < 0.0 ? 1 : 0;
        }
        keys.push_back("orientation: front " + rnd_hex(front) + ", back " + rnd_hex(back) + " (" +
                       std::to_string(back_faces) + " of " + std::to_string(num_faces) +
                       " drawn faces)");
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


    if (FixedFit) {
        // The bounding sphere stands in for the projected extents, so the model
        // keeps one size and one place while the camera moves.
        const double r = rPrep.mRadius;
        if (perspective) {
            min_x = min_y = -1.0;
            max_x = max_y = 1.0;
            min_z = 1.0 / (fit_distance + r);
            max_z = 1.0 / (fit_distance - r);
        } else {
            const double cx = rnd_dot(rPrep.mCentre.data(), cam.mU);
            const double cy = rnd_dot(rPrep.mCentre.data(), cam.mV);
            const double cz = rnd_dot(rPrep.mCentre.data(), cam.mW);
            min_x = cx - r;
            max_x = cx + r;
            min_y = cy - r;
            max_y = cy + r;
            min_z = cz - r;
            max_z = cz + r;
        }
        any = true;
    }

    // Fit the drawn geometry to the frame (leaving room for a colour bar).
    const int ss = rOpt.mSupersample;
    const double width = static_cast<double>(rOpt.mWidth) * ss;
    const double height = static_cast<double>(rOpt.mHeight) * ss;
    const double aspect = rOpt.mPixelAspect;
    const bool colorbar = rOpt.mColorbar && continuous;
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

    // Back faces seen through a cut-away are the inside: they take the tint.
    auto faces_away = [&](std::size_t f) {
        const double* nrm = &f_normal[3 * f];
        double toward[3] = {cam.mW[0], cam.mW[1], cam.mW[2]};
        if (perspective) {
            double c[3] = {0, 0, 0};
            const std::int64_t cnt = f_start[f + 1] - f_start[f];
            for (std::int64_t k = f_start[f]; k < f_start[f + 1]; ++k)
                for (std::size_t a = 0; a < 3; ++a)
                    c[a] += xyz[3 * static_cast<std::size_t>(s_nodes[static_cast<std::size_t>(k)]) + a];
            for (std::size_t a = 0; a < 3; ++a)
                toward[a] = eye[a] - (cnt > 0 ? c[a] / static_cast<double>(cnt) : 0.0);
        }
        return nrm[0] * toward[0] + nrm[1] * toward[1] + nrm[2] * toward[2] < 0.0;
    };

    // Triangles.
    std::size_t num_tris = 0;
    for (std::size_t f = 0; f < num_faces; ++f) {
        const std::int64_t n = f_start[f + 1] - f_start[f];
        if (n >= 3)
            num_tris += static_cast<std::size_t>(n - 2);
    }
    scene.mTris.reserve(num_tris);
    std::vector<std::array<int, 3>> tri_positions;
    for (std::size_t f = 0; f < num_faces; ++f) {
        const std::int64_t* ring = s_nodes.data() + f_start[f];
        const std::int64_t n = f_start[f + 1] - f_start[f];
        if (n < 3)
            continue;
        RndColor color = coloring.Map(f_values[f], rOpt.mFillColor);
        if (clip_on && !orientation && faces_away(f))
            color = rOpt.mCutawayTint;
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
        // The diagonal of a quad is chosen on the face's own (unsplit) ring, so
        // smooth shading's copies of a point never change it.
        rnd_face_triangles(xyz, f_nodes.data() + f_start[f], n, tri_positions);
        for (const std::array<int, 3>& t : tri_positions) {
            detail::RasterTri tri;
            for (std::size_t k = 0; k < 3; ++k) {
                tri.mV[k] = ring[t[k]];
                tri.mIntensity[k] = corner_i(ring[t[k]]);
            }
            tri.mColor = color;
            tri.mId = f_ids[f];
            scene.mTris.push_back(tri);
        }
    }

    // Lines: the input's line cells, then the edge overlays and extra layers.
    if (clip_on) {
        for (const RndClipped::Line& cut : clipped.mLines) {
            detail::RasterLine line;
            line.mA = cut.mA;
            line.mB = cut.mB;
            line.mColor = coloring.MapOther(l_values[cut.mSource], rOpt.mLineColor);
            line.mId = rPrep.mLines.mIds[cut.mSource];
            scene.mLines.push_back(line);
        }
        scene.mLines.insert(scene.mLines.end(), clipped.mEdgeLines.begin(),
                            clipped.mEdgeLines.end());
    } else {
        for (std::size_t i = 0; i < rPrep.mLines.mIds.size(); ++i) {
            detail::RasterLine line;
            line.mA = rPrep.mLines.mNodes[static_cast<std::size_t>(rPrep.mLines.mStart[i])];
            line.mB = rPrep.mLines.mNodes[static_cast<std::size_t>(rPrep.mLines.mStart[i]) + 1];
            line.mColor = coloring.MapOther(l_values[i], rOpt.mLineColor);
            line.mId = rPrep.mLines.mIds[i];
            scene.mLines.push_back(line);
        }
        scene.mLines.insert(scene.mLines.end(), rPrep.mEdgeLines.begin(), rPrep.mEdgeLines.end());
    }
    if (rPrep.mDeferCategoryEdges)
        rnd_category_edge_lines(rOpt, categorical, f_start, f_nodes, f_values, scene.mLines);
    for (const Extra& e : extras)
        scene.mLines.push_back({e.mA, e.mB, e.mColor, -1});
    auto add_point = [&](std::size_t i) {
        detail::RasterPoint point;
        point.mV = rPrep.mPoints.mNodes[static_cast<std::size_t>(rPrep.mPoints.mStart[i])];
        point.mColor = coloring.MapOther(p_values[i], rOpt.mLineColor);
        point.mId = rPrep.mPoints.mIds[i];
        scene.mPoints.push_back(point);
    };
    if (clip_on) {
        for (std::size_t i : clipped.mPoints)
            add_point(i);
    } else {
        for (std::size_t i = 0; i < rPrep.mPoints.mIds.size(); ++i)
            add_point(i);
    }

    const detail::RasterTarget target =
        detail::rasterize(scene, rOpt.mWidth * ss, rOpt.mHeight * ss, rOpt.mBackground);
    Frame frame;
    frame.mWidth = rOpt.mWidth;
    frame.mHeight = rOpt.mHeight;
    detail::raster_downsample(target, ss, frame.mRgba, frame.mCellIds);
    if (rOpt.mTheme == RenderTheme::Synthwave)
        rnd_theme_background(frame, rOpt);

    // Overlays and notes, at the output resolution.
    const double min_side = static_cast<double>(std::min(rOpt.mWidth, rOpt.mHeight));
    const std::int64_t margin =
        std::max<std::int64_t>(1, static_cast<std::int64_t>(min_side * 0.03));
    std::int64_t bar_w = 0;
    if (continuous) {
        frame.mColored = true;
        frame.mVMin = coloring.mVMin;
        frame.mVMax = coloring.mVMax;
        std::string text = label;
        if (component.has_value() && rOpt.mReduce.empty())
            text += "[" + std::to_string(*component) + "]";
        text +=
            ": " + rnd_num(coloring.mVMin) + " .. " + rnd_num(coloring.mVMax) + " (" + rOpt.mCmap;
        if (rOpt.mScale == RenderScale::Log)
            text += ", log";
        else if (rOpt.mScale == RenderScale::Symlog)
            text += ", symlog " + rnd_num(rOpt.mScaleThreshold);
        text += ")";
        frame.mNotes.push_back(text);
        std::size_t non_finite = 0;
        for (const std::vector<double>* p_list : {&f_values, &l_values, &p_values})
            for (double v : *p_list)
                if (!std::isfinite(v) || (rOpt.mScale == RenderScale::Log && !(v > 0.0)))
                    ++non_finite;
        if (rOpt.mDiagnostic == RenderDiagnostic::EdgeLength)
            non_finite -= std::min(non_finite, l_values.size() + p_values.size());
        if (non_finite > 0)
            frame.mNotes.push_back(
                std::to_string(non_finite) + " drawn " + (non_finite == 1 ? "value" : "values") +
                " with no finite value on this scale, in " + rnd_hex(rOpt.mNanColor));
    }
    if (colorbar) {
        bar_w = std::max<std::int64_t>(2, rOpt.mWidth / 40);
        const std::int64_t x1 = rOpt.mWidth - 1 - margin;
        const std::int64_t x0 = x1 - bar_w + 1;
        const std::int64_t y0 = static_cast<std::int64_t>(rOpt.mHeight * 0.2);
        const std::int64_t y1 =
            std::max<std::int64_t>(y0 + 1, static_cast<std::int64_t>(rOpt.mHeight * 0.8));
        const double s_lo =
            detail::render_scale_forward(rOpt.mScale, rOpt.mScaleThreshold, coloring.mVMin);
        const double s_hi =
            detail::render_scale_forward(rOpt.mScale, rOpt.mScaleThreshold, coloring.mVMax);
        for (std::int64_t y = y0; y <= y1; ++y) {
            const double t = static_cast<double>(y1 - y) / static_cast<double>(y1 - y0);
            // The bar is linear in the scaled value, so its colour at height t is
            // the colour of the scaled value s_lo + t * (s_hi - s_lo).
            const detail::Rgb c = detail::colormap_lookup(
                coloring.mpTable, std::isfinite(s_lo) && s_hi > s_lo ? t : 0.5);
            for (std::int64_t x = x0; x <= x1; ++x)
                rnd_plot(frame, x, y, {c.mR, c.mG, c.mB, 255});
        }
        // Ticks beside the bar, and the same values as a line of text.
        const std::vector<double> ticks = detail::render_legend_ticks(
            rOpt.mScale, rOpt.mScaleThreshold, coloring.mVMin, coloring.mVMax, 5);
        std::string tick_text;
        const RndColor tick_color = {110, 110, 110, 255};
        for (double v : ticks) {
            const double s = detail::render_scale_forward(rOpt.mScale, rOpt.mScaleThreshold, v);
            if (!std::isfinite(s) || !(s_hi > s_lo))
                continue;
            const double t = (s - s_lo) / (s_hi - s_lo);
            const double y = static_cast<double>(y1) - t * static_cast<double>(y1 - y0);
            rnd_line(frame, static_cast<double>(x0 - 3), y, static_cast<double>(x0 - 1), y,
                     tick_color);
            tick_text += (tick_text.empty() ? "" : ", ") + rnd_num(v);
        }
        if (!tick_text.empty())
            frame.mNotes.push_back("ticks: " + tick_text);
    }
    for (const std::string& key : keys)
        frame.mNotes.push_back(key);
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
    rnd_theme_post(frame, rOpt);
    return frame;
}

}  // namespace

RenderScene prepare_render(const Mesh& rMesh, const RenderOptions& rOpt) {
    RenderScene scene;
    scene.mpData = rnd_prepare(rMesh, rOpt);
    return scene;
}

bool render_scene_bounds(const RenderScene& rScene, std::array<double, 3>& rCentre,
                         double& rRadius) {
    if (!rScene.mpData)
        return false;
    const RndPrepared& prep = *static_cast<const RndPrepared*>(rScene.mpData.get());
    rCentre = prep.mCentre;
    rRadius = prep.mRadius;
    return true;
}

Frame render_scene(const RenderScene& rScene, const RenderOptions& rOpt, bool FixedFit) {
    if (!rScene.mpData)
        throw std::invalid_argument(std::string(kRndPrefix) + "render_scene needs a prepared scene");
    return rnd_draw(*static_cast<const RndPrepared*>(rScene.mpData.get()), rOpt, FixedFit);
}

Frame render(const Mesh& rMesh, const RenderOptions& rOpt) {
    const std::shared_ptr<const RndPrepared> prepared = rnd_prepare(rMesh, rOpt);
    return rnd_draw(*prepared, rOpt, false);
}

}  // namespace meshioplusplus
