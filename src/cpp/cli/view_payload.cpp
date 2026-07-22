// SPDX-License-Identifier: MIT
/// @file view_payload.cpp
/// @brief Implementation of the mesh -> Polyscope mapping. Twin of
///        `src/python/meshioplusplus/_viewer.py` -- KEEP THEM IN SYNC.

#include "view_payload.hpp"
#include "polyscope_view.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus::cli {

namespace {

/// Volume cell types Polyscope holds directly. Everything else 3D is
/// simplexified into tetrahedra first.
bool view_is_native_volume_type(const std::string& rType) {
    return rType == "tetra" || rType == "hexahedron";
}

/// Components above which per-component scalars stop being useful.
constexpr std::size_t kMaxComponentQuantities = 16;

/// Topological dimension of a cell type, including the ragged prefixes the
/// shared table does not key directly.
int view_dimension_of(const std::string& rType) {
    if (rType.rfind("polyhedron", 0) == 0)
        return 3;
    if (rType.rfind("polygon", 0) == 0)
        return 2;
    const auto& table = topological_dimension();
    auto it = table.find(rType);
    return it == table.end() ? -1 : it->second;
}

int view_max_dimension(const Mesh& rMesh) {
    int best = 0;
    for (const auto& block : rMesh.CellRange())
        best = std::max(best, view_dimension_of(block.Type()));
    return best;
}

/// `mesh.points` as `(P, 3)`, padded or truncated with a note.
std::vector<std::array<double, 3>> view_vertices_of(const Mesh& rMesh,
                                                    std::vector<std::string>& rNotes) {
    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    std::vector<std::array<double, 3>> out(n, {0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t c = 0; c < std::min<std::size_t>(dim, 3); ++c)
            out[i][c] = detail::read_double(points, i * dim + c);

    if (dim < 3)
        rNotes.push_back("points are " + std::to_string(dim) +
                         "D; padded to 3D with zeros for display");
    else if (dim > 3)
        rNotes.push_back("points are " + std::to_string(dim) +
                         "D; only the first 3 coordinates are displayed");
    return out;
}

/**
 * @brief The blocks of one topological dimension, plus their global cell ids.
 *
 * `rCellSource[i]` is the index of emitted primitive `i` in the *global* cell
 * numbering across every block. The counter advances over **skipped** blocks
 * too -- that is the whole point. A gmsh file whose `triangle` block follows a
 * boundary-marker `line` block is what makes it necessary: a naive
 * per-dimension concatenation offsets every cell quantity by the length of the
 * skipped block, and the result still renders.
 */
std::vector<std::size_t> view_dimension_rows(const Mesh& rMesh, int Dim,
                                             std::vector<std::size_t>& rBlockIndices) {
    std::vector<std::size_t> source;
    std::size_t global = 0;
    std::size_t index = 0;
    for (const auto& block : rMesh.CellRange()) {
        const std::size_t n = block.NumCells();
        if (view_dimension_of(block.Type()) == Dim) {
            rBlockIndices.push_back(index);
            for (std::size_t c = 0; c < n; ++c)
                source.push_back(global + c);
        }
        global += n;
        ++index;
    }
    return source;
}

/// One cell's connectivity, as vertex indices.
std::vector<std::uint32_t> view_row(const Mesh::CellView& rBlock, std::size_t Cell) {
    std::vector<std::uint32_t> out;
    if (rBlock.IsRagged()) {
        const std::size_t n = rBlock.RowSize(Cell);
        out.reserve(n);
        for (std::size_t k = 0; k < n; ++k)
            out.push_back(static_cast<std::uint32_t>(rBlock.Row(Cell)[k]));
        return out;
    }
    const std::size_t npc = rBlock.NodesPerCell();
    out.reserve(npc);
    const NDArray& conn = rBlock.Conn();
    for (std::size_t k = 0; k < npc; ++k)
        out.push_back(static_cast<std::uint32_t>(
            detail::read_double(conn, Cell * npc + k)));
    return out;
}

/**
 * @brief `cell_data[name]` flattened to one value per *global* cell.
 *
 * Returns false (with a note) when the array cannot be interpreted: meshio++
 * stores one sub-array per block, and a malformed or inconsistent set has no
 * single flat meaning. Refusing beats guessing -- a wrong answer here is a
 * plausible-looking but wrongly-coloured picture.
 */
bool view_flatten_cell_data(const Mesh& rMesh, const std::string& rName,
                            std::vector<double>& rOut, std::size_t& rComponents,
                            std::vector<std::string>& rNotes) {
    const std::size_t blocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(rName) != blocks) {
        rNotes.push_back("cell_data '" + rName + "' has " +
                         std::to_string(rMesh.CellDataNumBlocks(rName)) + " block(s) for " +
                         std::to_string(blocks) + " cell block(s); dropped");
        return false;
    }

    rOut.clear();
    rComponents = 0;
    std::size_t b = 0;
    for (const auto& block : rMesh.CellRange()) {
        const NDArray& a = rMesh.CellData(rName, b);
        const std::size_t n = block.NumCells();
        const std::size_t comps = n == 0 ? 1 : a.Size() / std::max<std::size_t>(n, 1);
        if (b == 0)
            rComponents = comps;
        if (comps != rComponents) {
            rNotes.push_back("cell_data '" + rName +
                             "' blocks disagree on component count; dropped");
            return false;
        }
        if (a.Size() != n * comps) {
            rNotes.push_back("cell_data '" + rName + "' block length does not match the " +
                             std::to_string(n) + "-cell '" + block.Type() +
                             "' block; dropped");
            return false;
        }
        for (std::size_t i = 0; i < a.Size(); ++i)
            rOut.push_back(detail::read_double(a, i));
        ++b;
    }
    if (rComponents == 0)
        rComponents = 1;
    return true;
}

/// Whether a 3-wide array should render as RGB rather than as arrows.
///
/// Gated on the **name**, never inferred from the value range: a normalized
/// displacement field or a unit normal shifted into [0,1] would otherwise be
/// silently drawn as colour with no way for the user to tell.
bool view_looks_like_color(const std::string& rName, const std::vector<double>& rValues) {
    auto lower = rName;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    auto ends_with = [&lower](const std::string& rSuffix) {
        return lower.size() >= rSuffix.size() &&
               lower.compare(lower.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
    };
    const bool named = lower == "color" || lower == "colour" || lower == "colors" ||
                       lower == "colours" || lower == "rgb" || ends_with("_color") ||
                       ends_with("_colour") || ends_with("_colors") ||
                       ends_with("_colours") || ends_with("_rgb") || ends_with(":rgb");
    if (!named)
        return false;

    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (double v : rValues) {
        if (!std::isfinite(v))
            continue;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }
    if (lo > hi)
        return false;
    return lo >= -1e-9 && hi <= 1.0 + 1e-9;
}

/// Finite range of a value list, or false when every entry is non-finite.
bool view_finite_range(const std::vector<double>& rValues, double& rMin, double& rMax,
                       bool& rHasNonFinite) {
    rMin = std::numeric_limits<double>::infinity();
    rMax = -std::numeric_limits<double>::infinity();
    rHasNonFinite = false;
    for (double v : rValues) {
        if (!std::isfinite(v)) {
            rHasNonFinite = true;
            continue;
        }
        rMin = std::min(rMin, v);
        rMax = std::max(rMax, v);
    }
    return rMin <= rMax;
}

/// Map one data array onto zero or more quantities.
void view_quantities_from_array(const std::string& rName, const std::vector<double>& rFlat,
                                std::size_t Components, QuantityOn On,
                                const std::string& rSource,
                                std::vector<Quantity>& rOut,
                                std::vector<std::string>& rNotes) {
    const std::size_t n = Components == 0 ? 0 : rFlat.size() / Components;

    auto component = [&](std::size_t c) {
        std::vector<double> out(n);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = rFlat[i * Components + c];
        return out;
    };
    auto magnitude = [&]() {
        std::vector<double> out(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t c = 0; c < Components; ++c) {
                const double v = rFlat[i * Components + c];
                sum += v * v;
            }
            out[i] = std::sqrt(sum);
        }
        return out;
    };
    auto push_scalar = [&](const std::string& rQName, std::vector<double> Values) {
        double lo = 0.0;
        double hi = 0.0;
        bool nonFinite = false;
        if (!view_finite_range(Values, lo, hi, nonFinite)) {
            rNotes.push_back(rSource + " '" + rQName + "' has no finite values; dropped");
            return;
        }
        Quantity q;
        q.mName = rQName;
        q.mKind = QuantityKind::Scalar;
        q.mOn = On;
        q.mValues = std::move(Values);
        // Only pin the range when non-finite entries are present; otherwise
        // leave Polyscope to autoscale as it normally would.
        q.mHasRange = nonFinite;
        q.mMin = lo;
        q.mMax = hi;
        rOut.push_back(std::move(q));
    };
    auto push_wide = [&](const std::string& rQName, QuantityKind Kind,
                         std::vector<double> Values) {
        Quantity q;
        q.mName = rQName;
        q.mKind = Kind;
        q.mOn = On;
        q.mValues = std::move(Values);
        rOut.push_back(std::move(q));
    };

    if (Components == 1) {
        push_scalar(rName, component(0));
        return;
    }
    if (Components == 2) {
        rNotes.push_back(rSource + " '" + rName + "' is 2-wide; zero-padded to a 3D vector");
        std::vector<double> padded(n * 3, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            padded[i * 3] = rFlat[i * 2];
            padded[i * 3 + 1] = rFlat[i * 2 + 1];
        }
        push_wide(rName, QuantityKind::Vector, std::move(padded));
        return;
    }
    if (Components == 3) {
        if (view_looks_like_color(rName, rFlat)) {
            std::vector<double> rgb = rFlat;
            for (double& v : rgb)
                v = std::min(1.0, std::max(0.0, v));
            push_wide(rName, QuantityKind::Color, std::move(rgb));
            return;
        }
        push_wide(rName, QuantityKind::Vector, rFlat);
        push_scalar(rName + ":magnitude", magnitude());
        return;
    }
    if (Components == 6 || Components == 9) {
        rNotes.push_back(rSource + " '" + rName + "' is a " + std::to_string(Components) +
                         "-component tensor; Polyscope has no tensor quantity, so "
                         "components and the Frobenius norm are shown");
    }
    if (Components > kMaxComponentQuantities) {
        rNotes.push_back(rSource + " '" + rName + "' has " + std::to_string(Components) +
                         " components; only its magnitude is shown");
        push_scalar(rName + ":magnitude", magnitude());
        return;
    }
    for (std::size_t c = 0; c < Components; ++c)
        push_scalar(rName + ":" + std::to_string(c), component(c));
    if (Components == 6 || Components == 9)
        push_scalar(rName + ":norm", magnitude());
}

/// Gather every usable `cell_data` array onto the emitted primitives.
void view_collect_quantities(const Mesh& rPointSource, const Mesh& rCellSource,
                             const std::vector<std::size_t>& rCellIds, QuantityOn PointOn,
                             QuantityOn CellOn, std::vector<Quantity>& rOut,
                             std::vector<std::string>& rNotes) {
    for (const auto& name : rPointSource.PointDataNames()) {
        const NDArray& a = rPointSource.PointData(name);
        const std::size_t n = rPointSource.NumPoints();
        const std::size_t comps = n == 0 ? 1 : a.Size() / std::max<std::size_t>(n, 1);
        std::vector<double> flat(a.Size());
        for (std::size_t i = 0; i < a.Size(); ++i)
            flat[i] = detail::read_double(a, i);
        view_quantities_from_array(name, flat, comps, PointOn, "point_data", rOut, rNotes);
    }

    std::set<std::string> pointNames;
    for (const auto& q : rOut)
        pointNames.insert(q.mName);

    for (const auto& name : rCellSource.CellDataNames()) {
        std::vector<double> flat;
        std::size_t comps = 1;
        if (!view_flatten_cell_data(rCellSource, name, flat, comps, rNotes))
            continue;

        // Provenance is plumbing, not something to colour by.
        if (name == "surface:parent_cell")
            continue;

        std::vector<double> gathered(rCellIds.size() * comps);
        for (std::size_t i = 0; i < rCellIds.size(); ++i)
            for (std::size_t c = 0; c < comps; ++c)
                gathered[i * comps + c] = flat[rCellIds[i] * comps + c];

        // Point and cell arrays share one Polyscope name space per structure,
        // so a name present on both sides would silently shadow.
        std::string qname = name;
        if (pointNames.count(name)) {
            qname = name + " (cells)";
            rNotes.push_back("'" + name +
                             "' exists as both point and cell data; the cell one is shown "
                             "as '" +
                             qname + "'");
        }
        view_quantities_from_array(qname, gathered, comps, CellOn, "cell_data", rOut,
                                   rNotes);
    }
}

/// Whether linearizing actually changed the mesh, so the note is only emitted
/// when something was really lost.
bool view_types_changed(const Mesh& rBefore, const Mesh& rAfter) {
    if (rBefore.NumPoints() != rAfter.NumPoints())
        return true;
    auto before = rBefore.CellRange().begin();
    auto after = rAfter.CellRange().begin();
    const auto beforeEnd = rBefore.CellRange().end();
    const auto afterEnd = rAfter.CellRange().end();
    for (; before != beforeEnd && after != afterEnd; ++before, ++after)
        if ((*before).Type() != (*after).Type())
            return true;
    return false;
}

ConvertCellsOptions view_convert_options(ConvertCellsMode Mode) {
    ConvertCellsOptions options;
    options.mMode = Mode;
    return options;
}

}  // namespace

std::size_t ViewPayload::NumPrimitives() const {
    if (!mFaces.empty())
        return mFaces.size();
    if (!mEdges.empty())
        return mEdges.size();
    if (!mTets.empty())
        return mTets.size();
    if (!mHexes.empty())
        return mHexes.size();
    return mMixedCells.size();
}

ViewKind view_kind_from_name(const std::string& rName) {
    if (rName == "auto")
        return ViewKind::Auto;
    if (rName == "surface")
        return ViewKind::Surface;
    if (rName == "volume")
        return ViewKind::Volume;
    if (rName == "curve")
        return ViewKind::Curve;
    if (rName == "points")
        return ViewKind::Points;
    throw std::invalid_argument("view: unknown kind '" + rName +
                                "' (expected auto, surface, volume, curve or points)");
}

namespace {

ViewPayload view_build_points(const Mesh& rMesh, std::vector<std::string> Notes) {
    ViewPayload payload;
    payload.mKind = ViewKind::Points;
    payload.mVertices = view_vertices_of(rMesh, Notes);
    view_collect_quantities(rMesh, rMesh, {}, QuantityOn::Vertices, QuantityOn::Vertices,
                            payload.mQuantities, Notes);
    payload.mNotes = std::move(Notes);
    return payload;
}

ViewPayload view_build_curve(const Mesh& rMesh, std::vector<std::string> Notes) {
    ViewPayload payload;
    payload.mKind = ViewKind::Curve;
    payload.mVertices = view_vertices_of(rMesh, Notes);

    std::vector<std::size_t> blocks;
    const std::vector<std::size_t> ids = view_dimension_rows(rMesh, 1, blocks);
    bool trimmed = false;
    for (std::size_t b : blocks) {
        const auto block = rMesh.Cells(b);
        for (std::size_t c = 0; c < block.NumCells(); ++c) {
            const auto row = view_row(block, c);
            if (row.size() > 2)
                trimmed = true;
            if (row.size() >= 2)
                payload.mEdges.push_back({row[0], row[1]});
        }
    }
    if (trimmed)
        Notes.push_back("higher-order line cells are drawn as their end segments");

    view_collect_quantities(rMesh, rMesh, ids, QuantityOn::Vertices, QuantityOn::Edges,
                            payload.mQuantities, Notes);
    payload.mNotes = std::move(Notes);
    return payload;
}

ViewPayload view_build_surface(const Mesh& rMesh, std::vector<std::string> Notes) {
    // Linearize rather than truncating columns by hand: convert_cells preserves
    // cell data 1:1, prunes the now-orphaned mid-side nodes, and remaps point
    // data -- all of which would otherwise have to be redone here.
    //
    // Run unconditionally. The linear-base table lives in convert_cells.cpp's
    // anonymous namespace, and duplicating it here to decide whether the call
    // is needed would create a *second* thing to keep in sync; convert_cells is
    // idempotent on types it does not apply to, so the only cost is one copy of
    // an already-linear mesh, which a viewer can afford.
    Mesh linearized =
        convert_cells(rMesh, view_convert_options(ConvertCellsMode::Linearize)).mMesh;
    if (view_types_changed(rMesh, linearized))
        Notes.push_back("higher-order surface cells are drawn as their linear bases");
    const Mesh* pSource = &linearized;

    ViewPayload payload;
    payload.mKind = ViewKind::Surface;
    payload.mVertices = view_vertices_of(*pSource, Notes);

    std::vector<std::size_t> blocks;
    const std::vector<std::size_t> ids = view_dimension_rows(*pSource, 2, blocks);
    for (std::size_t b : blocks) {
        const auto block = pSource->Cells(b);
        for (std::size_t c = 0; c < block.NumCells(); ++c)
            payload.mFaces.push_back(view_row(block, c));
    }

    view_collect_quantities(*pSource, *pSource, ids, QuantityOn::Vertices,
                            QuantityOn::Faces, payload.mQuantities, Notes);
    payload.mNotes = std::move(Notes);
    return payload;
}

ViewPayload view_build_surface_of_volume(const Mesh& rMesh, std::vector<std::string> Notes) {
    Mesh surface = extract_surface(rMesh, /*recordParentIds=*/true);
    gather_cell_data_onto_surface(rMesh, surface);
    if (!rMesh.CellDataNames().empty())
        Notes.push_back("cell data is sampled from each boundary face's owning volume cell");

    Mesh linear =
        convert_cells(surface, view_convert_options(ConvertCellsMode::Linearize)).mMesh;
    if (view_types_changed(surface, linear))
        Notes.push_back("the extracted boundary is drawn with linearized cells");
    surface = std::move(linear);

    ViewPayload payload;
    payload.mKind = ViewKind::Surface;
    payload.mVertices = view_vertices_of(surface, Notes);

    std::vector<std::size_t> blocks;
    const std::vector<std::size_t> ids = view_dimension_rows(surface, 2, blocks);
    for (std::size_t b : blocks) {
        const auto block = surface.Cells(b);
        for (std::size_t c = 0; c < block.NumCells(); ++c)
            payload.mFaces.push_back(view_row(block, c));
    }

    view_collect_quantities(surface, surface, ids, QuantityOn::Vertices, QuantityOn::Faces,
                            payload.mQuantities, Notes);
    payload.mNotes = std::move(Notes);
    return payload;
}

ViewPayload view_build_volume(const Mesh& rMesh, std::vector<std::string> Notes) {
    std::set<std::string> types;
    for (const auto& block : rMesh.CellRange())
        if (view_dimension_of(block.Type()) == 3)
            types.insert(block.Type());

    for (const auto& t : types)
        if (t.rfind("polyhedron", 0) == 0)
            throw std::invalid_argument(
                "view: cannot render polyhedron cells ('" + t +
                "') as a volume mesh; pass --kind surface to draw the boundary instead");

    const Mesh* pSource = &rMesh;
    Mesh simplexified;
    bool allNative = true;
    for (const auto& t : types)
        if (!view_is_native_volume_type(t))
            allNative = false;

    if (!allNative) {
        // Polyscope holds only tets and hexes, so wedges/pyramids/higher-order
        // must be decomposed. Any hexes present are collateral damage.
        std::string note =
            "volume cells are not renderable directly and were split into tetrahedra";
        if (types.count("hexahedron"))
            note += " -- hexahedra in the same mesh were split too";
        Notes.push_back(note);
        simplexified =
            convert_cells(rMesh, view_convert_options(ConvertCellsMode::Simplexify)).mMesh;
        pSource = &simplexified;
        types.clear();
        for (const auto& block : pSource->CellRange())
            if (view_dimension_of(block.Type()) == 3)
                types.insert(block.Type());
    }

    ViewPayload payload;
    payload.mKind = ViewKind::Volume;
    payload.mVertices = view_vertices_of(*pSource, Notes);

    std::vector<std::size_t> blocks;
    const std::vector<std::size_t> ids = view_dimension_rows(*pSource, 3, blocks);

    const bool onlyTets = types.size() == 1 && types.count("tetra") == 1;
    const bool onlyHexes = types.size() == 1 && types.count("hexahedron") == 1;

    for (std::size_t b : blocks) {
        const auto block = pSource->Cells(b);
        for (std::size_t c = 0; c < block.NumCells(); ++c) {
            const auto row = view_row(block, c);
            if (onlyTets) {
                payload.mTets.push_back({row[0], row[1], row[2], row[3]});
            } else if (onlyHexes) {
                std::array<std::uint32_t, 8> hex{};
                for (std::size_t k = 0; k < 8; ++k)
                    hex[k] = row[k];
                payload.mHexes.push_back(hex);
            } else {
                // Ours, in block order -- see the note on ViewPayload::mMixedCells.
                std::array<std::uint32_t, 8> cell{};
                cell.fill(kPadIndex);
                for (std::size_t k = 0; k < std::min<std::size_t>(row.size(), 8); ++k)
                    cell[k] = row[k];
                payload.mMixedCells.push_back(cell);
            }
        }
    }

    view_collect_quantities(*pSource, *pSource, ids, QuantityOn::Vertices,
                            QuantityOn::Cells, payload.mQuantities, Notes);
    payload.mNotes = std::move(Notes);
    return payload;
}

}  // namespace

ViewPayload build_view_payload(const Mesh& rMesh, ViewKind Kind) {
    std::vector<std::string> notes;
    const int dim = view_max_dimension(rMesh);

    if (Kind == ViewKind::Points || (Kind == ViewKind::Auto && dim == 0))
        return view_build_points(rMesh, std::move(notes));
    if (Kind == ViewKind::Curve || (Kind == ViewKind::Auto && dim == 1))
        return view_build_curve(rMesh, std::move(notes));
    if (Kind == ViewKind::Volume) {
        if (dim != 3)
            throw std::invalid_argument(
                "view: --kind volume needs 3D cells, but the mesh has none (highest cell "
                "dimension is " +
                std::to_string(dim) + "); use --kind surface");
        return view_build_volume(rMesh, std::move(notes));
    }
    if (Kind == ViewKind::Surface && dim == 3)
        return view_build_surface_of_volume(rMesh, std::move(notes));
    if (Kind == ViewKind::Auto && dim == 3)
        return view_build_volume(rMesh, std::move(notes));
    return view_build_surface(rMesh, std::move(notes));
}

}  // namespace meshioplusplus::cli

#ifndef MESHIOPLUSPLUS_HAS_POLYSCOPE

// Always defined, so `view`/`screenshot` exist in every build and fail with the
// flag that enables them rather than not existing at all -- the same contract
// partition_kahip_parts follows for KaHIP. Never a link error, never a silent
// no-op. The real implementations live in polyscope_view.cpp, which is the only
// translation unit that includes a Polyscope header.

namespace meshioplusplus::cli {

namespace {
[[noreturn]] void view_not_built(const char* pVerb) {
    throw std::runtime_error(
        std::string(pVerb) +
        ": this build has no viewer; rebuild with -DMESHIOPLUSPLUS_WITH_POLYSCOPE=ON "
        "(and `git submodule update --init --recursive`)");
}
}  // namespace

bool has_polyscope() noexcept { return false; }

void view_mesh(const Mesh&, ViewKind, const std::string&, const std::string&) {
    view_not_built("view");
}

void screenshot_mesh(const Mesh&, ViewKind, const std::string&, const std::string&,
                     const std::string&, int, int, bool) {
    view_not_built("screenshot");
}

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_HAS_POLYSCOPE
