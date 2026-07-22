#include "meshioplusplus/detail/face_color.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// "no point/cell data array named 'x' (available: a, b)" -- the writers'
// equivalent of data_unknown_key_message, kept local because this one has to
// name both locations at once.
std::string fc_unknown_message(const Mesh& rMesh, const std::string& rName) {
    std::string msg = "meshio++: no point_data or cell_data array named '" + rName + "'";
    std::string avail;
    for (const std::string& n : rMesh.PointDataNames())
        avail += (avail.empty() ? "" : ", ") + ("point:" + n);
    for (const std::string& n : rMesh.CellDataNames())
        avail += (avail.empty() ? "" : ", ") + ("cell:" + n);
    if (avail.empty())
        return msg + " (the mesh carries no data arrays)";
    return msg + " (available: " + avail + ")";
}

// Components = product of the trailing dimensions. A 1-D array is one scalar
// per entry; anything else spreads its remaining axes across the row.
std::size_t fc_num_components(const NDArray& rArray, std::size_t NumEntries) {
    if (NumEntries == 0)
        return 0;
    const std::size_t total = rArray.Size();
    return total / NumEntries;
}

// Reduce one row to a scalar: the requested component, or the magnitude.
//
// The magnitude sums left to right and takes one sqrt at the end -- the Python
// twin does the same, rather than np.linalg.norm, so the two round alike.
double fc_scalarize(const NDArray& rArray, std::size_t Row, std::size_t NumComponents,
                    const std::optional<int>& rComponent) {
    if (NumComponents == 0)
        return std::nan("");
    const std::size_t base = Row * NumComponents;
    if (rComponent.has_value()) {
        const int c = *rComponent;
        if (c < 0 || static_cast<std::size_t>(c) >= NumComponents)
            throw std::invalid_argument("meshio++: component " + std::to_string(c) +
                                        " is out of range for an array with " +
                                        std::to_string(NumComponents) + " component(s)");
        return read_double(rArray, base + static_cast<std::size_t>(c));
    }
    if (NumComponents == 1)
        return read_double(rArray, base);
    double sum = 0.0;
    for (std::size_t k = 0; k < NumComponents; ++k) {
        const double v = read_double(rArray, base + k);
        sum += v * v;
    }
    return std::sqrt(sum);
}

// Flatten a cell_data array to one scalar per global cell of `rMesh`, in the
// same block-major order mSourceCell and "surface:parent_cell" count in.
std::vector<double> fc_flatten_cell_data(const Mesh& rMesh, const std::string& rName,
                                         const std::optional<int>& rComponent) {
    const std::size_t num_blocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(rName) != num_blocks)
        throw std::invalid_argument("meshio++: cell_data array '" + rName + "' has " +
                                    std::to_string(rMesh.CellDataNumBlocks(rName)) +
                                    " block(s) but the mesh has " + std::to_string(num_blocks));
    std::vector<double> out;
    for (std::size_t b = 0; b < num_blocks; ++b) {
        const std::size_t n = rMesh.Cells(b).NumCells();
        const NDArray& arr = rMesh.CellData(rName, b);
        const std::size_t ncomp = fc_num_components(arr, n);
        for (std::size_t r = 0; r < n; ++r)
            out.push_back(fc_scalarize(arr, r, ncomp, rComponent));
    }
    return out;
}

// The drawn mesh's per-face parent cell in the source mesh, flattened in the
// drawn mesh's own block-major cell order -- so indexing it by a face's
// mSourceCell yields the source cell that owns it. Empty when the drawn mesh
// carries no provenance (it is the source mesh, undrilled).
std::vector<std::int64_t> fc_flatten_parent_ids(const Mesh& rDraw) {
    const std::string key = "surface:parent_cell";
    const std::size_t num_blocks = rDraw.NumCellBlocks();
    if (!rDraw.HasCellData(key) || rDraw.CellDataNumBlocks(key) != num_blocks)
        return {};
    std::vector<std::int64_t> out;
    for (std::size_t b = 0; b < num_blocks; ++b) {
        const std::size_t n = rDraw.Cells(b).NumCells();
        const NDArray& arr = rDraw.CellData(key, b);
        for (std::size_t r = 0; r < n; ++r)
            out.push_back(read_int(arr, r));
    }
    return out;
}

}  // namespace

std::optional<Rgb> FaceColors::Color(std::size_t Face) const {
    if (!mActive || Face >= mValues.size())
        return std::nullopt;
    const double v = mValues[Face];
    if (!std::isfinite(v))
        return std::nullopt;
    return colormap_lookup(mpTable, color_param(v, mVMin, mVMax));
}

double color_param(double v, double vMin, double vMax) {
    // An all-constant array gives vMin == vMax; the middle of the colormap is a
    // better answer than a division by zero, and matches the Python twin.
    if (!(vMax > vMin))
        return 0.5;
    const double t = (v - vMin) / (vMax - vMin);
    if (!(t > 0.0))
        return 0.0;
    if (t > 1.0)
        return 1.0;
    return t;
}

std::vector<ColorFace> color_faces_from_projection(const std::vector<ProjectedFace>& rFaces) {
    std::vector<ColorFace> out;
    out.reserve(rFaces.size());
    for (const ProjectedFace& f : rFaces)
        out.push_back(ColorFace{f.mNodes, f.mNumNodes, f.mSourceCell});
    return out;
}

std::vector<ColorFace> color_faces_flat(const Mesh& rMesh) {
    std::vector<ColorFace> out;
    // Mirrors both flat writers' loop: block-major over the drawable types,
    // then cell-major, with the global cell counter advancing over skipped
    // blocks too (the mSourceCell convention).
    std::int64_t global_cell_base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t block_base = global_cell_base;
        global_cell_base += static_cast<std::int64_t>(cb.NumCells());
        const std::string& type = cb.Type();
        if (type != "line" && type != "triangle" && type != "quad")
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t ncols = cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            ColorFace face;
            face.mNodes = {-1, -1, -1, -1};
            face.mNumNodes = static_cast<std::uint8_t>(ncols < 4 ? ncols : 4);
            for (std::size_t k = 0; k < face.mNumNodes; ++k)
                face.mNodes[k] = read_int(conn, r * ncols + k);
            face.mSourceCell = block_base + static_cast<std::int64_t>(r);
            out.push_back(face);
        }
    }
    return out;
}

FaceColors resolve_face_colors(const ColorSpec& rSpec, const Mesh& rSource, const Mesh& rDraw,
                               const std::vector<ColorFace>& rFaces) {
    FaceColors out;
    if (!rSpec.Active())
        return out;

    if (rSpec.mVMin.has_value() && rSpec.mVMax.has_value() && *rSpec.mVMin > *rSpec.mVMax)
        throw std::invalid_argument("meshio++: vmin must not exceed vmax");

    // Resolve the colormap first: an unknown name should fail before any work.
    out.mpTable = colormap_table(rSpec.mCmap);

    const std::string& name = rSpec.mColorBy;
    const bool is_point = rSource.HasPointData(name);
    if (!is_point && !rSource.HasCellData(name))
        throw std::invalid_argument(fc_unknown_message(rSource, name));

    out.mValues.reserve(rFaces.size());

    if (is_point) {
        // Point data lives on the drawn mesh: a boundary extraction carries
        // point_data through its point compaction, so the drawn mesh's node ids
        // index the already-subset array.
        if (!rDraw.HasPointData(name))
            throw std::invalid_argument("meshio++: point_data array '" + name +
                                        "' did not survive boundary extraction");
        const NDArray& arr = rDraw.PointData(name);
        const std::size_t ncomp = fc_num_components(arr, rDraw.NumPoints());
        for (const ColorFace& f : rFaces) {
            // Mean of the corner values, summed in node order then divided
            // once -- the same shape as the projection's depth accumulation,
            // and what the Python twin reproduces (never np.mean, whose
            // pairwise summation rounds differently).
            double sum = 0.0;
            for (std::uint8_t k = 0; k < f.mNumNodes; ++k)
                sum += fc_scalarize(arr, static_cast<std::size_t>(f.mNodes[k]), ncomp,
                                    rSpec.mComponent);
            out.mValues.push_back(f.mNumNodes == 0 ? std::nan("")
                                                   : sum / static_cast<double>(f.mNumNodes));
        }
    } else {
        const std::vector<double> flat = fc_flatten_cell_data(rSource, name, rSpec.mComponent);
        const std::vector<std::int64_t> parents = fc_flatten_parent_ids(rDraw);
        for (const ColorFace& f : rFaces) {
            std::int64_t cell = f.mSourceCell;
            if (!parents.empty()) {
                if (cell < 0 || static_cast<std::size_t>(cell) >= parents.size()) {
                    out.mValues.push_back(std::nan(""));
                    continue;
                }
                cell = parents[static_cast<std::size_t>(cell)];
            }
            if (cell < 0 || static_cast<std::size_t>(cell) >= flat.size())
                out.mValues.push_back(std::nan(""));
            else
                out.mValues.push_back(flat[static_cast<std::size_t>(cell)]);
        }
    }

    // Auto range over the DRAWN faces, not the whole array: the visible figure
    // then spans the whole colorbar. (ParaView's default is the whole array;
    // this is documented as differing.) Non-finite values are excluded, with
    // FiniteStats' seed-on-first-finite rule rather than +/-DBL_MAX sentinels.
    double lo = 0.0;
    double hi = 0.0;
    bool seen = false;
    for (double v : out.mValues) {
        if (!std::isfinite(v))
            continue;
        if (!seen) {
            lo = v;
            hi = v;
            seen = true;
        } else {
            if (v < lo)
                lo = v;
            if (v > hi)
                hi = v;
        }
    }
    out.mVMin = rSpec.mVMin.has_value() ? *rSpec.mVMin : lo;
    out.mVMax = rSpec.mVMax.has_value() ? *rSpec.mVMax : hi;
    if (out.mVMin > out.mVMax)
        throw std::invalid_argument("meshio++: vmin must not exceed vmax");
    out.mActive = true;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
