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

// System includes
#include <algorithm>
#include <cmath>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

NDArray data_owned_copy(const NDArray& rArray) {
    NDArray c = rArray;
    c.MakeOwned();
    return c;
}

Mesh clone_geometry(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(data_owned_copy(rMesh.Points()));
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                cells[c].resize(cb.NumFaces(c));
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    cells[c][f].assign(face.first, face.first + face.second);
                }
            }
            out.AddPolyhedronBlock(std::string(cb.Type()), std::move(cells));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                rows[c].assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            out.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
        } else {
            out.AddCellBlock(std::string(cb.Type()), data_owned_copy(cb.Conn()));
        }
    }
    // Named regions ride through verbatim. Every caller of this function leaves
    // the point and cell numbering untouched (the data operations, `transform`,
    // `smooth`, `interpolate`'s target clone, `attach_quality`), so there is
    // nothing to remap — see doc/regions.md's table of which operations remap
    // and which pass through.
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        out.AddRegion(rMesh.Region(i));
    // Property sets ride along too: these operations preserve the mesh's shape,
    // so dropping a deck's material data here would be new lossiness. They are
    // keyed by id, not by entity index, so there is nothing to remap.
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        out.AddPropertySet(rMesh.GetPropertySet(i));
    return out;
}

void accumulate_stats(const NDArray& rArray, std::size_t NumComponents,
                      std::vector<FiniteStats>& rStats) {
    if (rStats.size() < NumComponents)
        rStats.resize(NumComponents);
    const std::size_t total = rArray.Size();
    if (total == 0 || NumComponents == 0)
        return;
    const std::size_t nrows = total / NumComponents;

    const std::size_t grain = 4096;
    const std::size_t nchunks = (nrows + grain - 1) / grain;
    std::vector<std::vector<FiniteStats>> partial(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t ci) {
            std::vector<FiniteStats> local(NumComponents);
            const std::size_t begin = ci * grain;
            const std::size_t end = std::min(begin + grain, nrows);
            for (std::size_t r = begin; r < end; ++r)
                for (std::size_t k = 0; k < NumComponents; ++k)
                    local[k].Add(read_double(rArray, r * NumComponents + k));
            partial[ci] = std::move(local);
        },
        1);
    for (const std::vector<FiniteStats>& chunk : partial)
        for (std::size_t k = 0; k < NumComponents && k < chunk.size(); ++k)
            rStats[k].Merge(chunk[k]);
}

FiniteStats combine_components(const std::vector<FiniteStats>& rStats) {
    FiniteStats all;
    for (const FiniteStats& s : rStats)
        all.Merge(s);
    return all;
}

namespace {

/**
 * @brief Unsigned area of a corner polygon (triangle / quad) via the Newell
 * normal. Mirrors `stats_area` in `operations/stats.cpp`. Internal to
 * `cell_measure`, so not exposed in the header.
 * @param rCoords the corner coordinates.
 * @param Corners how many of them form the polygon.
 * @return the unsigned area, or 0 for fewer than 3 corners.
 */
double dataops_polygon_area(const std::vector<Vec3>& rCoords, int Corners) {
    if (Corners < 3)
        return 0.0;
    Vec3 s = {0, 0, 0};
    for (int i = 0; i < Corners; ++i)
        s = vec3_add(s, vec3_cross(rCoords[i], rCoords[(i + 1) % Corners]));
    return 0.5 * vec3_norm(s);
}

/**
 * @brief Signed volume of a 3D cell via the divergence theorem over its
 * outward-wound boundary faces. Mirrors `stats_signed_volume`. Internal to
 * `cell_measure`, so not exposed in the header.
 * @param rCoords the cell's corner coordinates.
 * @param Type the cell type, which supplies the face table.
 * @return the signed volume, or NaN for a type with no face table.
 */
double dataops_cell_signed_volume(const std::vector<Vec3>& rCoords, CellType Type) {
    const std::vector<CellFaceDef>& faces = cell_faces(Type);
    if (faces.empty())
        return std::nan("");
    double vol6 = 0.0;
    for (const CellFaceDef& f : faces) {
        const Vec3 a = rCoords[f.mNodes[0]];
        for (int i = 1; i + 1 < f.mNumCorners; ++i)
            vol6 += triple_product(a, rCoords[f.mNodes[i]], rCoords[f.mNodes[i + 1]]);
    }
    return vol6 / 6.0;
}

}  // namespace

double cell_measure(const NDArray& rPoints, std::size_t PointDim, const Mesh::CellView& rCell,
                    std::size_t Index) {
    if (rCell.IsRagged() || rCell.IsPolyhedron())
        return std::nan("");
    const CellType ct = cell_type_from_name(rCell.Type());
    const int corners = cell_corner_count(ct);
    if (corners <= 0)
        return std::nan("");
    std::vector<Vec3> coords;
    read_corner_coords(rPoints, PointDim, rCell.Conn(), Index * rCell.NodesPerCell(),
                       static_cast<std::size_t>(corners), coords);
    const int dim = cell_type_dimension(ct);
    if (dim == 3) {
        const double v = dataops_cell_signed_volume(coords, ct);
        return std::isnan(v) ? v : std::fabs(v);
    }
    if (dim == 2)
        return dataops_polygon_area(coords, corners);
    if (dim == 1 && corners >= 2)
        return vec3_norm(vec3_sub(coords[1], coords[0]));
    return std::nan("");
}

}  // namespace detail
}  // namespace meshioplusplus
