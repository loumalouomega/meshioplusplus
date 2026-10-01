// SPDX-License-Identifier: MIT
#include "vtk_xml_pieces.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "region_field_data.hpp"
#include "text_cursor.hpp"
#include "vtk_xml_read.hpp"

namespace meshioplusplus::detail {
namespace {

std::size_t vxp_product(std::size_t a, std::size_t b) {
    if (b && a > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) / b)
        throw ReadError("VTK: structured dimensions overflow");
    return a * b;
}

std::size_t vxp_sum(std::size_t a, std::size_t b) {
    if (b > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) - a)
        throw ReadError("VTK: piece counts overflow");
    return a + b;
}

template <class T, std::size_t N>
std::array<T, N> vxp_attribute(pugi::xml_node node, const char* name,
                               std::array<T, N> fallback, bool required = false) {
    const auto attr = node.attribute(name);
    if (!attr) {
        if (required)
            throw ReadError(std::string("VTK: missing ") + name);
        return fallback;
    }
    TextStream in(attr.as_string());
    for (auto& value : fallback)
        if (!(in >> value))
            throw ReadError(std::string("VTK: malformed ") + name);
    std::string extra;
    if (in >> extra)
        throw ReadError(std::string("VTK: malformed ") + name);
    return fallback;
}

std::array<std::int64_t, 3> vxp_dimensions(const std::array<std::int64_t, 6>& rExtent) {
    std::array<std::int64_t, 3> dims;
    for (std::size_t k = 0; k < 3; ++k) {
        const auto lo = rExtent[2 * k], hi = rExtent[2 * k + 1];
        if (hi < lo || (lo < 0 && hi > std::numeric_limits<std::int64_t>::max() + lo - 1) ||
            hi - lo == std::numeric_limits<std::int64_t>::max())
            throw ReadError("VTK: inverted or overflowing extent");
        dims[k] = hi - lo + 1;
    }
    return dims;
}

struct VxpPiece {
    pugi::xml_node mNode;
    std::array<std::int64_t, 6> mExtent{};
    std::array<std::int64_t, 3> mDims{};
    std::size_t mNumPoints = 0;
    std::size_t mNumCells = 0;
};

struct VxpHeader {
    pugi::xml_node mGrid;
    VtuContext mCtx;
    std::vector<VxpPiece> mPieces;
    std::array<double, 3> mOrigin{};
    std::array<double, 3> mSpacing{{1, 1, 1}};
    bool mPoly = false;
    bool mImage = false;
};

VxpHeader vxp_header(const VtuSource& rSource, const char* pType, const char* pFormat) {
    const auto root = rSource.mDoc.child("VTKFile");
    if (!root || std::string(root.attribute("type").as_string()) != pType)
        throw ReadError(std::string("VTK: expected ") + pType);
    VxpHeader h;
    h.mGrid = root.child(pType);
    if (!h.mGrid)
        throw ReadError(std::string("VTK: missing ") + pType);
    h.mCtx = vtk_xml_read_context(rSource, pFormat);
    h.mPoly = std::string(pType) == "PolyData";
    h.mImage = std::string(pType) == "ImageData";
    std::array<std::int64_t, 6> whole{};
    if (!h.mPoly) {
        whole = vxp_attribute<std::int64_t, 6>(h.mGrid, "WholeExtent", {}, true);
        vxp_dimensions(whole);
    }
    if (h.mImage) {
        h.mOrigin = vxp_attribute<double, 3>(h.mGrid, "Origin", {});
        h.mSpacing = vxp_attribute<double, 3>(h.mGrid, "Spacing", {{1, 1, 1}});
        const std::array<double, 9> identity{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
        if (vxp_attribute<double, 9>(h.mGrid, "Direction", identity) != identity)
            throw ReadError("VTI with a non-identity Direction is not supported");
    }
    const bool many = h.mGrid.child("Piece").next_sibling("Piece");
    for (const auto node : h.mGrid.children("Piece")) {
        VxpPiece piece;
        piece.mNode = node;
        if (h.mPoly) {
            const std::string count = node.attribute("NumberOfPoints").as_string();
            if (count.empty() || count.find_first_not_of("0123456789") != std::string::npos)
                throw ReadError("VTP: invalid NumberOfPoints");
            try {
                piece.mNumPoints = vxp_sum(0, std::stoull(count));
            } catch (const std::exception&) {
                throw ReadError("VTP: invalid NumberOfPoints");
            }
        } else {
            piece.mExtent = vxp_attribute<std::int64_t, 6>(node, "Extent", whole, many);
            piece.mDims = vxp_dimensions(piece.mExtent);
            for (std::size_t k = 0; k < 3; ++k)
                if (piece.mExtent[2 * k] < whole[2 * k] ||
                    piece.mExtent[2 * k + 1] > whole[2 * k + 1])
                    throw ReadError("VTK: Piece Extent is outside WholeExtent");
            piece.mNumPoints = vxp_product(vxp_product(piece.mDims[0], piece.mDims[1]),
                                           piece.mDims[2]);
            piece.mNumCells = vxp_product(vxp_product(piece.mDims[0] - 1, piece.mDims[1] - 1),
                                          piece.mDims[2] - 1);
        }
        h.mPieces.push_back(piece);
    }
    if (h.mPieces.empty())
        throw ReadError("No Piece found");
    return h;
}

NDArray vxp_array(pugi::xml_node da, const VtuContext& rCtx, std::size_t rows,
                  int default_components = 1) {
    if (!da)
        throw ReadError("VTK: missing DataArray");
    int nc = 0;
    NDArray arr = vtu_read_data_array(da, rCtx, nc);
    if (!da.attribute("NumberOfComponents"))
        nc = default_components;
    if (arr.Size() != vxp_product(rows, nc))
        throw ReadError("VTK: DataArray length differs from piece count");
    if (nc > 1)
        arr.Reshape({rows, static_cast<std::size_t>(nc)});
    return arr;
}

bool vxp_same_shape(const NDArray& rA, const NDArray& rB) {
    return rA.Dtype() == rB.Dtype() && rA.Shape().size() == rB.Shape().size() &&
           std::equal(rA.Shape().begin() + 1, rA.Shape().end(), rB.Shape().begin() + 1);
}

NDArray vxp_concat(const std::vector<NDArray>& rParts) {
    auto shape = rParts.front().Shape();
    shape[0] = 0;
    for (const auto& arr : rParts) {
        if (!vxp_same_shape(rParts.front(), arr))
            throw ReadError("VTK: pieces disagree on array dtype or components");
        shape[0] = vxp_sum(shape[0], arr.Shape()[0]);
    }
    NDArray out = NDArray::Uninit(rParts.front().Dtype(), shape);
    std::size_t offset = 0;
    for (const auto& arr : rParts) {
        if (arr.Nbytes())
            std::memcpy(out.Data() + offset, arr.Data(), arr.Nbytes());
        offset += arr.Nbytes();
    }
    return out;
}

void vxp_poly_cells(const VxpPiece& rPiece, const VtuContext& rCtx,
                     std::vector<std::int64_t>& rConn, std::vector<std::int64_t>& rOffsets,
                     std::vector<std::int64_t>& rTypes, bool Metadata) {
    const std::array<const char*, 4> tags{{"Verts", "Lines", "Polys", "Strips"}};
    for (std::size_t kind = 0; kind < tags.size(); ++kind) {
        const auto section = rPiece.mNode.child(tags[kind]);
        if (!section)
            continue;
        std::vector<std::int64_t> conn, offsets;
        for (const auto da : section.children("DataArray")) {
            const std::string name = da.attribute("Name").as_string();
            if (name != "offsets" && (Metadata || name != "connectivity"))
                continue;
            int nc = 0;
            auto values = vtu_to_int64(vtu_read_data_array(da, rCtx, nc));
            if (name == "offsets")
                offsets = std::move(values);
            else
                conn = std::move(values);
        }
        std::int64_t prev = 0;
        const auto base = Metadata ? (rOffsets.empty() ? 0 : rOffsets.back())
                                   : static_cast<std::int64_t>(rConn.size());
        for (const auto end : offsets) {
            if (end <= prev)
                throw ReadError("VTP: non-increasing cell offsets");
            const auto size = end - prev;
            if (kind == 0 && size != 1)
                throw ReadError("poly-vertex VTP cells are not supported");
            if (kind == 1 && size != 2)
                throw ReadError("poly-line VTP cells are not supported");
            if (kind == 2 && size < 3)
                throw ReadError("VTP: polygon has fewer than three points");
            if (kind == 3)
                throw ReadError("triangle-strip VTP cells are not supported");
            rTypes.push_back(kind == 0 ? 1 : kind == 1 ? 3 : size == 3 ? 5 : size == 4 ? 9 : 7);
            rOffsets.push_back(vxp_sum(base, end));
            prev = end;
        }
        if (!Metadata && static_cast<std::size_t>(prev) != conn.size())
            throw ReadError("VTP: offsets do not span connectivity");
        for (const auto index : conn)
            if (index < 0 || static_cast<std::size_t>(index) >= rPiece.mNumPoints)
                throw ReadError("VTP: point index out of range");
        rConn.insert(rConn.end(), conn.begin(), conn.end());
    }
}

std::vector<std::pair<std::string, NDArray>> vxp_fields(pugi::xml_node holder,
                                                       const VtuContext& rCtx,
                                                       const ReadOptions& rOpts) {
    std::vector<std::pair<std::string, NDArray>> out;
    for (const auto da : holder.child("FieldData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!is_region_field_name(name) && (!rOpts.WantsAnyData() || !rOpts.WantsArray(name)))
            continue;
        try {
            dtype_from_vtu(da.attribute("type").as_string());
        } catch (const ReadError&) {
            log::warn("VTK: skipping non-numeric field array '{}'", name);
            continue;
        }
        int nc = 0;
        NDArray arr = vtu_read_data_array(da, rCtx, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
        out.emplace_back(name, std::move(arr));
    }
    return out;
}

using VxpData = std::map<std::string, std::vector<NDArray>>;

void vxp_data(const VxpPiece& rPiece, const VtuContext& rCtx, const ReadOptions& rOpts,
               const char* pSection, std::size_t count, VxpData& rOut) {
    if (!rOpts.WantsAnyData())
        return;
    for (const auto da : rPiece.mNode.child(pSection).children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (rOpts.WantsArray(name))
            rOut[name].push_back(vxp_array(da, rCtx, count));
    }
}

std::vector<std::string> vxp_names(const VxpHeader& rHeader, const char* pSection) {
    // Metadata follows the same cross-piece dtype/component compatibility rule.
    std::map<std::string, std::vector<std::pair<std::string, int>>> names;
    for (const auto& piece : rHeader.mPieces)
        for (const auto da : piece.mNode.child(pSection).children("DataArray"))
            names[da.attribute("Name").as_string()].emplace_back(
                da.attribute("type").as_string(), da.attribute("NumberOfComponents").as_int(1));
    std::vector<std::string> out;
    for (const auto& [name, parts] : names)
        if (parts.size() == rHeader.mPieces.size() &&
            std::all_of(parts.begin(), parts.end(), [&](const auto& part) { return part == parts[0]; }))
            out.push_back(name);
    return out;
}

}  // namespace

void vtk_structured_cells(const std::array<std::int64_t, 3>& rDims,
                           std::vector<std::int64_t>& rConn,
                           std::vector<std::int64_t>& rOffsets,
                           std::vector<std::int64_t>& rTypes, bool VolumeOnly) {
    std::vector<std::size_t> axes;
    std::array<std::int64_t, 3> cells;
    for (std::size_t k = 0; k < 3; ++k) {
        if (rDims[k] < 1)
            throw ReadError("VTK: dimensions must be positive");
        if (rDims[k] > 1)
            axes.push_back(k);
        cells[k] = std::max<std::int64_t>(rDims[k] - 1, 1);
    }
    if (VolumeOnly && axes.size() != 3)
        return;
    const std::array<std::int64_t, 3> stride{{1, rDims[0],
        static_cast<std::int64_t>(vxp_product(rDims[0], rDims[1]))}};
    const auto count = vxp_product(vxp_product(cells[0], cells[1]), cells[2]);
    const std::size_t width = std::size_t{1} << axes.size();
    vxp_product(count, width);
    for (std::int64_t k = 0; k < cells[2]; ++k)
        for (std::int64_t j = 0; j < cells[1]; ++j)
            for (std::int64_t i = 0; i < cells[0]; ++i) {
                const auto base = i + j * stride[1] + k * stride[2];
                const auto a = axes.empty() ? 0 : stride[axes[0]];
                const auto b = axes.size() < 2 ? 0 : stride[axes[1]];
                const auto c = axes.size() < 3 ? 0 : stride[axes[2]];
                const std::array<std::int64_t, 8> nodes{{base, base + a, base + a + b, base + b,
                    base + c, base + a + c, base + a + b + c, base + b + c}};
                rConn.insert(rConn.end(), nodes.begin(), nodes.begin() + width);
                rOffsets.push_back(rConn.size());
                rTypes.push_back(axes.empty() ? 1 : axes.size() == 1 ? 3 : axes.size() == 2 ? 9 : 12);
            }
}

Mesh vtk_xml_read_pieces(const std::string& rPath, const ReadOptions& rOpts,
                         const char* pType, const char* pFormat) {
    VtuSource source;
    vtu_load(rPath, pugi::parse_default, source, pType, pFormat);
    auto h = vxp_header(source, pType, pFormat);
    std::vector<NDArray> points;
    VxpData point_data, cell_data;
    std::vector<std::int64_t> conn, offsets, types;
    auto fields = vxp_fields(h.mGrid, h.mCtx, rOpts);
    std::map<std::string, NDArray> passthrough;
    std::map<std::string, std::vector<NDArray>> piece_regions;
    std::size_t point_base = 0, cell_base = 0;
    for (auto& piece : h.mPieces) {
        NDArray pts;
        if (h.mPoly || std::string(pType) == "StructuredGrid") {
            const auto da = piece.mNode.child("Points").child("DataArray");
            if (!da && piece.mNumPoints == 0)
                pts = NDArray(DType::Float64, {0, 3});
            else
                pts = vxp_array(da, h.mCtx, piece.mNumPoints, 3);
            if (pts.Shape().size() != 2 || pts.Shape()[1] != 3)
                throw ReadError("VTK: points must have three components");
        } else {
            std::array<std::vector<double>, 3> coordinates;
            if (!h.mImage) {
                std::vector<pugi::xml_node> arrays;
                for (const auto da : piece.mNode.child("Coordinates").children("DataArray"))
                    arrays.push_back(da);
                if (arrays.size() != 3)
                    throw ReadError("VTR Coordinates must have three axis DataArrays");
                const std::array<const char*, 3> names{{"x_coordinates", "y_coordinates", "z_coordinates"}};
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    auto selected = arrays[axis];
                    for (const auto da : arrays)
                        if (std::string(da.attribute("Name").as_string()) == names[axis])
                            selected = da;
                    const auto arr = vxp_array(selected, h.mCtx, piece.mDims[axis]);
                    if (arr.Shape().size() != 1)
                        throw ReadError("VTR: axis coordinates must be scalar");
                    for (std::size_t n = 0; n < arr.Size(); ++n)
                        coordinates[axis].push_back(read_double(arr, n));
                }
            }
            pts = NDArray::Uninit(DType::Float64, {piece.mNumPoints, 3});
            auto dst = pts.As<double>();
            std::size_t row = 0;
            for (std::int64_t k = 0; k < piece.mDims[2]; ++k)
                for (std::int64_t j = 0; j < piece.mDims[1]; ++j)
                    for (std::int64_t i = 0; i < piece.mDims[0]; ++i, ++row) {
                        const std::array<std::int64_t, 3> index{{i, j, k}};
                        for (std::size_t axis = 0; axis < 3; ++axis)
                            dst[row * 3 + axis] = h.mImage
                                ? (h.mOrigin[axis] + static_cast<double>(piece.mExtent[2 * axis]) * h.mSpacing[axis]) +
                                    static_cast<double>(index[axis]) * h.mSpacing[axis]
                                : coordinates[axis][index[axis]];
                    }
        }
        std::vector<std::int64_t> pc, po, pt;
        if (h.mPoly)
            vxp_poly_cells(piece, h.mCtx, pc, po, pt, false);
        else
            vtk_structured_cells(piece.mDims, pc, po, pt, true);
        piece.mNumCells = pt.size();
        vxp_data(piece, h.mCtx, rOpts, "PointData", piece.mNumPoints, point_data);
        vxp_data(piece, h.mCtx, rOpts, "CellData", pt.size(), cell_data);

        // Decode piece-local regions with the shared validator before shifting.
        // Dataset-level regions already index the assembled, piece-major file.
        auto local_fields = vxp_fields(piece.mNode, h.mCtx, rOpts);
        if (!local_fields.empty()) {
            Mesh local;
            local.AssignPoints(NDArray(pts));
            std::unordered_map<std::string, NDArray> no_data;
            std::vector<std::int64_t> file_map;
            reconstruct_cells(pc.data(), po, pt, no_data, nullptr, {}, local, &file_map);
            regions_from_field_arrays(local, local_fields, &file_map, pFormat);
            std::vector<std::int64_t> inverse(file_map.size());
            for (std::size_t n = 0; n < file_map.size(); ++n)
                inverse[file_map[n]] = n;
            auto encoded = regions_to_field_arrays(local, &inverse);
            for (auto& [name, arr] : encoded) {
                if (name.rfind("region-meta:", 0) == 0) {
                    fields.emplace_back(name, std::move(arr));
                    continue;
                }
                const bool point = name.rfind("region:point:", 0) == 0;
                const std::size_t step = name.rfind("region:side:", 0) == 0 ? 2 : 1;
                for (std::size_t n = 0; n < arr.Size(); n += step)
                    arr.As<std::int64_t>()[n] += point ? point_base : cell_base;
                piece_regions[name].push_back(std::move(arr));
            }
            for (const auto& name : local.FieldDataNames()) {
                if (is_region_field_name(name))
                    passthrough.insert_or_assign(name, NDArray(local.FieldData(name)));
                else
                    fields.emplace_back(name, NDArray(local.FieldData(name)));
            }
        }
        const auto conn_base = conn.size();
        for (const auto index : pc)
            conn.push_back(index + point_base);
        for (const auto end : po)
            offsets.push_back(vxp_sum(conn_base, end));
        types.insert(types.end(), pt.begin(), pt.end());
        point_base = vxp_sum(point_base, piece.mNumPoints);
        cell_base = vxp_sum(cell_base, pt.size());
        points.push_back(std::move(pts));
    }
    Mesh mesh;
    mesh.AssignPoints(vxp_concat(points));
    std::unordered_map<std::string, NDArray> raw_data;
    auto assemble = [&](VxpData& data, bool point) {
        for (auto& [name, parts] : data) {
            if (parts.size() != h.mPieces.size() ||
                !std::all_of(parts.begin(), parts.end(), [&](const auto& arr) {
                    return vxp_same_shape(parts.front(), arr);
                })) {
                log::warn("{}: data '{}' is missing from, or differs between, pieces; dropped", pFormat, name);
                continue;
            }
            auto arr = vxp_concat(parts);
            if (point)
                mesh.AddPointData(name, std::move(arr));
            else
                raw_data.emplace(name, std::move(arr));
        }
    };
    assemble(point_data, true);
    assemble(cell_data, false);
    check_vtk_cell_arrays(conn.size(), offsets, types, raw_data);
    std::vector<std::int64_t> file_map;
    reconstruct_cells(conn.data(), offsets, types, raw_data, nullptr, {}, mesh, &file_map);
    for (auto& [name, parts] : piece_regions)
        fields.emplace_back(name, vxp_concat(parts));
    regions_from_field_arrays(mesh, fields, &file_map, pFormat);
    for (auto& [name, arr] : passthrough)
        mesh.AddFieldData(name, std::move(arr));
    return mesh;
}

MeshMetadata vtk_xml_pieces_metadata(const std::string& rPath, const char* pType,
                                     const char* pFormat) {
    VtuSource source;
    vtu_load(rPath, pugi::parse_minimal, source, pType, pFormat);
    const auto h = vxp_header(source, pType, pFormat);
    MeshMetadata meta;
    meta.mPointDim = 3;
    meta.mHasBBox = h.mImage;
    std::vector<std::int64_t> conn, offsets, types;
    for (const auto& piece : h.mPieces) {
        meta.mNumPoints = vxp_sum(meta.mNumPoints, piece.mNumPoints);
        if (h.mPoly)
            vxp_poly_cells(piece, h.mCtx, conn, offsets, types, true);
        else if (piece.mNumCells) {
            // No geometry or array payload is decoded for structured metadata.
            CellBlockInfo info;
            info.mType = "hexahedron";
            info.mNumCells = piece.mNumCells;
            info.mNodesPerCell = 8;
            if (meta.mCellBlocks.empty())
                meta.mCellBlocks.push_back(info);
            else
                meta.mCellBlocks[0].mNumCells = vxp_sum(meta.mCellBlocks[0].mNumCells, piece.mNumCells);
        }
        if (h.mImage)
            for (std::size_t k = 0; k < 3; ++k) {
                const double a = h.mOrigin[k] + static_cast<double>(piece.mExtent[2 * k]) * h.mSpacing[k];
                const double b = a + static_cast<double>(piece.mDims[k] - 1) * h.mSpacing[k];
                const bool first = &piece == &h.mPieces.front();
                meta.mBBoxMin[k] = first ? std::min(a, b) : std::min(meta.mBBoxMin[k], std::min(a, b));
                meta.mBBoxMax[k] = first ? std::max(a, b) : std::max(meta.mBBoxMax[k], std::max(a, b));
            }
    }
    if (h.mPoly)
        meta.mCellBlocks = summarize_cells(offsets, types);
    meta.mPointDataNames = vxp_names(h, "PointData");
    meta.mCellDataNames = vxp_names(h, "CellData");
    std::vector<pugi::xml_node> holders{h.mGrid};
    for (const auto& piece : h.mPieces)
        holders.push_back(piece.mNode);
    std::map<std::string, std::size_t> region_counts;
    for (const auto holder : holders)
        for (const auto da : holder.child("FieldData").children("DataArray")) {
            const std::string name = da.attribute("Name").as_string();
            if (is_region_field_name(name)) {
                if (name.rfind("region:", 0) == 0)
                    region_counts[name] = vxp_sum(region_counts[name], da.attribute("NumberOfTuples").as_ullong());
            } else {
                try {
                    dtype_from_vtu(da.attribute("type").as_string());
                    meta.mFieldDataNames.push_back(name);
                } catch (const ReadError&) {
                }
            }
        }
    std::vector<std::pair<std::string, std::size_t>> regions(region_counts.begin(), region_counts.end());
    meta.mRegions = region_summaries_from_field_names(regions);
    auto& names = meta.mFieldDataNames;
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return meta;
}

}  // namespace meshioplusplus::detail
