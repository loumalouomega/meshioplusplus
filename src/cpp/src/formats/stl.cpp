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
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/chunked_parse.hpp"
#include "../detail/exact_weld.hpp"
#include "../detail/open_source.hpp"
#include "../detail/row_writer.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// First-occurrence de-duplication of 3-component rows, on their bytes.
// Returns per-row unique index; appends unique rows (raw bytes) to
// `rOutPoints`. The rows of an STL file sit at any alignment, so each
// coordinate's bytes are copied into its key word rather than read as a float.
std::vector<std::int64_t> dedup(const unsigned char* pRows, std::size_t nrows, std::size_t isz,
                                std::vector<unsigned char>& rOutPoints) {
    const std::size_t rowbytes = 3 * isz;
    std::vector<detail::ExactWeldKey> keys(nrows);
    parallel_for(nrows, [&](std::size_t i) {
        detail::ExactWeldKey& key = keys[i];
        key.mWidth = 3;
        for (std::size_t j = 0; j < 3; ++j)
            std::memcpy(&key.mWords[j], pRows + i * rowbytes + j * isz, isz);
    });
    detail::ExactWeld weld = detail::exact_weld_keys(keys);
    const std::size_t base = rOutPoints.size();
    rOutPoints.resize(base + weld.NumIds() * rowbytes);
    parallel_for(weld.NumIds(), [&](std::size_t id) {
        std::memcpy(rOutPoints.data() + base + id * rowbytes, pRows + weld.mFirstRow[id] * rowbytes,
                    rowbytes);
    });
    return std::move(weld.mIdOfRow);
}

Mesh build_mesh(std::vector<unsigned char>& rVertBytes, DType dt,
                std::vector<unsigned char>* pNormalBytes,
                const std::vector<std::uint16_t>* pAttributes = nullptr) {
    std::size_t isz = dtype_size(dt);
    std::size_t nverts = rVertBytes.size() / (3 * isz);

    std::vector<unsigned char> point_bytes;
    std::vector<std::int64_t> idx = dedup(rVertBytes.data(), nverts, isz, point_bytes);

    Mesh mesh;
    std::size_t num_unique = point_bytes.size() / (3 * isz);
    NDArray pts(dt, {num_unique, 3});
    if (!point_bytes.empty())
        std::memcpy(pts.Data(), point_bytes.data(), point_bytes.size());
    mesh.AssignPoints(std::move(pts));

    std::size_t ntri = nverts / 3;
    // An empty STL has no cells (match the Python reader, which returns no
    // cell blocks rather than an empty triangle block).
    if (ntri == 0)
        return mesh;

    NDArray cells(DType::Int64, {ntri, 3});
    std::int64_t* cp = cells.As<std::int64_t>();
    for (std::size_t i = 0; i < ntri * 3; ++i)
        cp[i] = idx[i];
    mesh.AddCellBlock("triangle", std::move(cells));

    if (pNormalBytes && !pNormalBytes->empty()) {
        NDArray nrm(DType::Float64, {ntri, 3});
        std::memcpy(nrm.Data(), pNormalBytes->data(), pNormalBytes->size());
        mesh.AppendCellData("facet_normals", std::move(nrm));
    }
    // The attribute word of a binary facet, kept only when some facet sets it.
    if (pAttributes && std::any_of(pAttributes->begin(), pAttributes->end(),
                                   [](std::uint16_t A) { return A != 0; })) {
        NDArray attr(DType::UInt16, {ntri});
        std::memcpy(attr.Data(), pAttributes->data(), ntri * sizeof(std::uint16_t));
        mesh.AppendCellData("stl:attribute", std::move(attr));
    }
    return mesh;
}

bool stl_starts_with(std::string_view S, std::string_view Prefix) {
    return S.substr(0, Prefix.size()) == Prefix;
}

bool stl_is_comment_line(std::string_view S) {
    return stl_starts_with(S, "solid") || stl_starts_with(S, "outer loop") ||
           stl_starts_with(S, "endloop") || stl_starts_with(S, "endfacet") ||
           stl_starts_with(S, "endsolid");
}

std::string_view stl_lstrip(std::string_view S) {
    std::size_t b = 0;
    while (b < S.size() && detail::text_is_blank(S[b]))
        ++b;
    return S.substr(b);
}

// The last 3 numbers of every non-comment line of `lines`, appended to rOut.
void stl_ascii_rows(const std::vector<std::string_view>& rLines, std::vector<double>& rOut) {
    // A line holds at most one row, so the line count bounds the rows.
    rOut.reserve(rOut.size() + rLines.size() * 3);
    std::vector<std::string_view> tok;
    for (const std::string_view line : rLines) {
        const std::string_view s = stl_lstrip(line);
        if (s.empty() || stl_is_comment_line(s))
            continue;
        detail::split_blanks(s, tok);
        if (tok.size() < 3)
            continue;
        for (std::size_t j = tok.size() - 3; j < tok.size(); ++j)
            rOut.push_back(detail::parse_double_prefix(tok[j]));
    }
}

struct StlSolid {
    std::string mName;
    std::size_t mFirst = 0;  // the first facet it holds
};

// Several `solid` records of an ASCII file become cell regions named after
// them (`solid_<k>` when unnamed, `_<k>` appended to a repeated name); one
// solid leaves the mesh as it was. `Full` is the whole file, `NumFacets` the
// facets the parse found, so a solid's range ends at the next one's start.
void stl_add_solid_regions(Mesh& rMesh, std::string_view Full, std::size_t NumFacets) {
    std::vector<StlSolid> solids;
    std::size_t data_lines = 0;
    std::size_t pos = 0;
    while (pos < Full.size()) {
        std::size_t eol = Full.find('\n', pos);
        if (eol == std::string_view::npos)
            eol = Full.size();
        const std::string_view line = stl_lstrip(Full.substr(pos, eol - pos));
        pos = eol + 1;
        if (line.empty())
            continue;
        if (stl_starts_with(line, "solid")) {
            std::string_view name = line.substr(5);
            while (!name.empty() && detail::text_is_blank(name.back()))
                name.remove_suffix(1);
            name = stl_lstrip(name);
            solids.push_back({std::string(name), data_lines / 4});
        } else if (!stl_is_comment_line(line)) {
            ++data_lines;
        }
    }
    if (solids.size() < 2)
        return;
    std::vector<std::string> used;
    for (std::size_t k = 0; k < solids.size(); ++k) {
        const std::size_t last = k + 1 < solids.size() ? solids[k + 1].mFirst : NumFacets;
        std::string name = solids[k].mName.empty() ? "solid_" + std::to_string(k) : solids[k].mName;
        while (std::find(used.begin(), used.end(), name) != used.end())
            name += "_" + std::to_string(k);
        used.push_back(name);
        const std::size_t count = last > solids[k].mFirst ? last - solids[k].mFirst : 0;
        NDArray entries = NDArray::Uninit(DType::Int64, {count});
        for (std::size_t q = 0; q < count; ++q)
            entries.As<std::int64_t>()[q] = static_cast<std::int64_t>(solids[k].mFirst + q);
        rMesh.AddRegion(Region(name, RegionKind::Cell, std::move(entries)));
    }
}

// `Text` is the ASCII body: the whole file, or the file after its header line.
// The views end where the file ends, so the numbers are parsed with the
// bounded prefix parse, never over a view that lacks a terminator.
Mesh read_ascii(std::string_view Text, std::string_view Full) {
    // Collect the last 3 numbers of every non-comment line; rows 0,4,8,... are
    // facet normals, the rest are vertices. A large file is read in chunks of
    // whole lines in parallel (roadmap 5.7.2), each chunk by the serial rule
    // into its own rows, joined in chunk order: the rows a serial read gives.
    std::vector<double> data;
    if (Text.size() >= detail::chunked_lines_min_bytes && detail::parallel_workers() > 1) {
        const std::vector<std::size_t> cuts = detail::chunk_lines(Text, 0, Text.size());
        std::vector<std::vector<double>> parts(cuts.size() - 1);
        parallel_for(
            parts.size(),
            [&](std::size_t c) {
                stl_ascii_rows(detail::split_lines(Text.substr(cuts[c], cuts[c + 1] - cuts[c])),
                               parts[c]);
            },
            1);
        std::size_t total = 0;
        for (const auto& rPart : parts)
            total += rPart.size();
        data.reserve(total);
        for (const auto& rPart : parts)
            data.insert(data.end(), rPart.begin(), rPart.end());
    } else {
        stl_ascii_rows(detail::split_lines(Text), data);
    }
    std::size_t nrows = data.size() / 3;
    if (nrows % 4 != 0)
        throw ReadError("Malformed ascii STL");

    // Facet f is rows 4f (its normal) and 4f+1..4f+3 (its corners).
    const std::size_t nfacets = nrows / 4;
    constexpr std::size_t kRow = 3 * sizeof(double);
    std::vector<unsigned char> verts(nfacets * 3 * kRow), normals(nfacets * kRow);
    parallel_for(nfacets, [&](std::size_t f) {
        const double* facet = data.data() + f * 12;
        std::memcpy(normals.data() + f * kRow, facet, kRow);
        std::memcpy(verts.data() + f * 3 * kRow, facet + 3, 3 * kRow);
    });
    Mesh mesh = build_mesh(verts, DType::Float64, &normals);
    if (mesh.NumCellBlocks() > 0)
        stl_add_solid_regions(mesh, Full, nfacets);
    return mesh;
}

// `pTriangles` points at the first of `num_tri` 50-byte records, whose
// presence the caller has checked against the file size.
Mesh read_binary(const char* pTriangles, std::uint32_t num_tri) {
    std::vector<unsigned char> verts;
    verts.reserve(static_cast<std::size_t>(num_tri) * 9 * sizeof(float));
    std::vector<std::uint16_t> attributes(num_tri);
    for (std::uint32_t i = 0; i < num_tri; ++i) {
        // bytes [12, 48) are the 9 float32 vertex coords (host is little-endian).
        const char* tri = pTriangles + static_cast<std::size_t>(i) * 50;
        verts.insert(verts.end(), reinterpret_cast<const unsigned char*>(tri + 12),
                     reinterpret_cast<const unsigned char*>(tri + 48));
        std::memcpy(&attributes[i], tri + 48, 2);
    }
    return build_mesh(verts, DType::Float32, nullptr, &attributes);
}

}  // namespace

Mesh read_stl(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "Could not open file: " + rPath);
    const std::string_view text = source.View();

    if (text.size() < 80)
        return read_ascii(text, text);

    // A file of 80..83 bytes has no complete triangle count: it is ASCII.
    if (text.size() >= 84) {
        std::uint32_t num_tri = 0;
        std::memcpy(&num_tri, text.data() + 80, 4);  // little-endian host
        if (84 + std::uint64_t(num_tri) * 50 == text.size())
            return read_binary(text.data() + 84, num_tri);
    }

    // Fall back to ascii: skip the first line.
    const std::size_t eol = text.find('\n');
    return read_ascii(eol == std::string_view::npos ? std::string_view() : text.substr(eol + 1),
                      text);
}

namespace {

void gather_triangles(const Mesh& rMesh, std::vector<std::array<double, 9>>& rTris,
                      std::vector<std::array<double, 3>>& rNormals,
                      std::vector<std::uint16_t>& rAttributes,
                      std::vector<std::int64_t>& rGlobal) {
    const bool have_attr = rMesh.HasCellData("stl:attribute");
    const bool have_normals = rMesh.HasCellData("facet_normals");
    const std::size_t normal_blocks = have_normals ? rMesh.CellDataNumBlocks("facet_normals") : 0;
    const NDArray& points = rMesh.Points();
    std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 3;
    const detail::DoubleView point_values(points);

    std::size_t block = 0;
    std::int64_t global_base = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle") {
            global_base += static_cast<std::int64_t>(cb.NumCells());
            ++block;
            continue;
        }
        std::size_t nc = cb.NumCells();
        const std::int64_t cell_base = global_base;
        global_base += static_cast<std::int64_t>(nc);
        const detail::Int64View conn(cb.Conn());
        std::optional<detail::DoubleView> attr;
        if (have_attr && block < rMesh.CellDataNumBlocks("stl:attribute"))
            attr.emplace(rMesh.CellData("stl:attribute", block));
        std::optional<detail::DoubleView> nrm;
        if (have_normals && block < normal_blocks)
            nrm.emplace(rMesh.CellData("facet_normals", block));
        for (std::size_t r = 0; r < nc; ++r) {
            std::array<double, 9> tri{};
            double v[3][3];
            for (int k = 0; k < 3; ++k) {
                std::int64_t pi = conn[r * 3 + k];
                for (int c = 0; c < 3; ++c)
                    v[k][c] = (std::size_t(c) < dim) ? point_values[pi * dim + c] : 0.0;
                tri[k * 3 + 0] = v[k][0];
                tri[k * 3 + 1] = v[k][1];
                tri[k * 3 + 2] = v[k][2];
            }
            rTris.push_back(tri);
            rGlobal.push_back(cell_base + static_cast<std::int64_t>(r));
            rAttributes.push_back(attr ? static_cast<std::uint16_t>((*attr)[r]) : 0);

            std::array<double, 3> n{};
            if (nrm) {
                for (int c = 0; c < 3; ++c)
                    n[c] = (*nrm)[r * 3 + c];
            } else {
                double a[3] = {v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]};
                double b[3] = {v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]};
                n[0] = a[1] * b[2] - a[2] * b[1];
                n[1] = a[2] * b[0] - a[0] * b[2];
                n[2] = a[0] * b[1] - a[1] * b[0];
                double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                if (len > 0) {
                    n[0] /= len;
                    n[1] /= len;
                    n[2] /= len;
                }
            }
            rNormals.push_back(n);
        }
        ++block;
    }
}

// Extract the linearized boundary skin of the volume cells and triangulate
// its quads (quad -> (0,1,2), (0,2,3)) into a single triangle block.
Mesh stl_skin_triangles(const Mesh& rMesh) {
    const Mesh skin = extract_skin(rMesh, /*linearize=*/true);
    std::vector<std::int64_t> tri;
    for (const auto cb : skin.CellRange()) {
        const detail::Int64View conn(cb.Conn());
        const std::size_t nc = cb.NumCells();
        if (cb.Type() == "triangle") {
            for (std::size_t r = 0; r < nc * 3; ++r)
                tri.push_back(conn[r]);
        } else {  // quad
            for (std::size_t r = 0; r < nc; ++r) {
                const std::int64_t a = conn[r * 4 + 0];
                const std::int64_t b = conn[r * 4 + 1];
                const std::int64_t c = conn[r * 4 + 2];
                const std::int64_t d = conn[r * 4 + 3];
                tri.insert(tri.end(), {a, b, c, a, c, d});
            }
        }
    }
    Mesh out;
    NDArray pts = skin.Points();
    out.AssignPoints(std::move(pts));
    NDArray block = NDArray::Uninit(DType::Int64, {tri.size() / 3, 3});
    std::memcpy(block.Data(), tri.data(), tri.size() * sizeof(std::int64_t));
    out.AddCellBlock("triangle", std::move(block));
    return out;
}

// Number of non-volume cell blocks (the ones the skin path drops).
std::size_t stl_num_surface_blocks(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        if (cell_type_dimension(cell_type_from_name(cb.Type())) != 3)
            ++n;
    return n;
}

}  // namespace

void write_stl(const std::string& rPath, const Mesh& rMesh, bool binary, bool skin) {
    if (skin && has_skinnable_cells(rMesh)) {
        const std::size_t dropped = stl_num_surface_blocks(rMesh);
        if (dropped > 0) {
            log::warn(
                "STL: writing the extracted skin of the volume cells; {} pre-existing "
                "non-volume cell block(s) dropped (pass skin=false for the legacy behavior).",
                dropped);
            detail::provenance_note("cells-dropped",
                                    std::to_string(dropped) +
                                        " non-volume cell block(s) dropped in favour of the "
                                        "extracted skin of the volume cells");
        }
        write_stl(rPath, stl_skin_triangles(rMesh), binary, /*skin=*/false);
        return;
    }

    std::vector<std::array<double, 9>> tris;
    std::vector<std::array<double, 3>> normals;
    std::vector<std::uint16_t> attributes;
    std::vector<std::int64_t> global_ids;
    gather_triangles(rMesh, tris, normals, attributes, global_ids);

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    if (binary) {
        char header[80];
        std::memset(header, '\0', 80);
        const std::string msg = detail::provenance_lines(detail::SlotTier::Bounded)[0];
        std::memcpy(header, msg.data(), std::min(msg.size(), sizeof(header)));
        os.write(header, 80);
        std::uint32_t n = static_cast<std::uint32_t>(tris.size());
        os.write(reinterpret_cast<const char*>(&n), 4);
        for (std::size_t i = 0; i < tris.size(); ++i) {
            float buf[12];
            for (int c = 0; c < 3; ++c)
                buf[c] = static_cast<float>(normals[i][c]);
            for (int c = 0; c < 9; ++c)
                buf[3 + c] = static_cast<float>(tris[i][c]);
            os.write(reinterpret_cast<const char*>(buf), 48);
            std::uint16_t attr = attributes[i];
            os.write(reinterpret_cast<const char*>(&attr), 2);
        }
    } else {
        // One facet's seven lines, formatted in parallel chunks, one number at
        // a time (`CNumber`: `to_chars` for `%.17g`, roadmap 5.7.1).
        const detail::CNumber number;
        const auto wr3 = [&](std::string& rBuf, const char* pPrefix, const double* p) {
            rBuf += pPrefix;
            for (int c = 0; c < 3; ++c) {
                rBuf.push_back(' ');
                number.Append(rBuf, "%.17g", p[c]);
            }
            rBuf.push_back('\n');
        };
        // One solid of the facets `pOrder` lists (every facet when null).
        const auto write_solid = [&](const std::string& rName,
                                     const std::vector<std::size_t>* pOrder) {
            os << (rName.empty() ? "solid\n" : "solid " + rName + "\n");
            detail::write_row_chunks(
                os, pOrder ? pOrder->size() : tris.size(),
                [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                    for (std::size_t q = First; q < Last; ++q) {
                        const std::size_t i = pOrder ? (*pOrder)[q] : q;
                        wr3(rBuf, "facet normal", normals[i].data());
                        rBuf += " outer loop\n";
                        wr3(rBuf, "  vertex", &tris[i][0]);
                        wr3(rBuf, "  vertex", &tris[i][3]);
                        wr3(rBuf, "  vertex", &tris[i][6]);
                        rBuf += " endloop\nendfacet\n";
                    }
                });
            os << (rName.empty() ? "endsolid\n" : "endsolid " + rName + "\n");
        };
        // Two or more cell regions: one solid each, the facets no region holds last.
        std::vector<std::size_t> region_ids;
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
            if (rMesh.Region(r).mKind == RegionKind::Cell)
                region_ids.push_back(r);
        if (region_ids.size() < 2) {
            write_solid("", nullptr);
        } else {
            std::unordered_map<std::int64_t, std::size_t> position;
            for (std::size_t i = 0; i < global_ids.size(); ++i)
                position.emplace(global_ids[i], i);
            std::vector<char> held(tris.size(), 0);
            for (const std::size_t r : region_ids) {
                const auto& region = rMesh.Region(r);
                const detail::Int64View entries(region.mEntries);
                std::vector<std::size_t> order;
                for (std::size_t q = 0; q < region.NumEntries(); ++q) {
                    const auto it = position.find(entries[q]);
                    if (it != position.end()) {
                        order.push_back(it->second);
                        held[it->second] = 1;
                    }
                }
                std::sort(order.begin(), order.end());
                write_solid(region.mName, &order);
            }
            std::vector<std::size_t> rest;
            for (std::size_t i = 0; i < tris.size(); ++i)
                if (!held[i])
                    rest.push_back(i);
            if (!rest.empty())
                write_solid("", &rest);
        }
    }
}

}  // namespace meshioplusplus
