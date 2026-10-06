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
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/open_source.hpp"
#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// Split "<stem>.node" / "<stem>.ele" into the two sibling paths.
std::pair<std::string, std::string> node_ele_paths(const std::string& rPath, bool& rOk) {
    std::size_t dot = rPath.find_last_of('.');
    rOk = false;
    if (dot == std::string::npos)
        return {"", ""};
    std::string suffix = rPath.substr(dot);
    std::string stem = rPath.substr(0, dot);
    if (suffix == ".node" || suffix == ".ele") {
        rOk = true;
        return {stem + ".node", stem + ".ele"};
    }
    return {"", ""};
}

// First non-comment, non-blank line is the header; remaining non-comment
// tokens (whitespace-separated, across lines) are the data stream. The tokens
// are views into the text they were split from, which the caller keeps alive.
struct Parsed {
    std::vector<std::string_view> mHeader;
    std::vector<std::string_view> mData;
};

Parsed parse_text(std::string_view Text, const std::string& rPath) {
    Parsed p;
    bool have_header = false;
    std::vector<std::string_view> tokens;
    detail::TextCursor cursor(Text);
    while (!cursor.AtEnd()) {
        const std::string_view line = cursor.Line();
        // trim leading whitespace
        std::size_t s = 0;
        while (s < line.size() && detail::text_is_blank(line[s]))
            ++s;
        if (s >= line.size() || line[s] == '#')
            continue;
        if (!have_header) {
            detail::split_blanks(line, p.mHeader);
            have_header = true;
        } else {
            detail::split_blanks(line, tokens);
            p.mData.insert(p.mData.end(), tokens.begin(), tokens.end());
        }
    }
    if (!have_header)
        throw ReadError("TetGen: missing header line in " + rPath);
    return p;
}

}  // namespace

Mesh read_tetgen(const std::string& rPath) {
    bool ok = false;
    auto paths = node_ele_paths(rPath, ok);
    if (!ok)
        throw ReadError("TetGen: expected a .node or .ele file");
    const std::string& node_path = paths.first;
    const std::string& ele_path = paths.second;

    Mesh mesh;

    // ---- nodes ----
    const detail::FileSource node_source =
        detail::open_source(node_path, "Could not open file: " + node_path);
    const Parsed nf = parse_text(node_source.View(), node_path);
    if (nf.mHeader.size() < 4)
        throw ReadError("TetGen: malformed .node header");
    std::int64_t npoints = detail::strtoll_token(nf.mHeader[0]);
    int dim = static_cast<int>(detail::strtoll_token(nf.mHeader[1]));
    int num_attrs = static_cast<int>(detail::strtoll_token(nf.mHeader[2]));
    int num_bmarkers = static_cast<int>(detail::strtoll_token(nf.mHeader[3]));
    if (dim != 3)
        throw ReadError("TetGen: need 3D points");

    const int ncol = 4 + num_attrs + num_bmarkers;
    if (static_cast<std::int64_t>(nf.mData.size()) != npoints * ncol)
        throw ReadError("TetGen: .node data size mismatch");

    auto at = [&](std::int64_t r, int c) -> double {
        return detail::parse_double_prefix(nf.mData[r * ncol + c]);
    };

    std::int64_t node_index_base = npoints > 0 ? static_cast<std::int64_t>(at(0, 0)) : 0;
    for (std::int64_t i = 0; i < npoints; ++i) {
        if (static_cast<std::int64_t>(at(i, 0)) != node_index_base + i)
            throw ReadError("TetGen: nodes not numbered consecutively");
    }

    NDArray pts(DType::Float64, {static_cast<std::size_t>(npoints), 3});
    double* pp = pts.As<double>();
    for (std::int64_t i = 0; i < npoints; ++i)
        for (int c = 0; c < 3; ++c)
            pp[i * 3 + c] = at(i, 1 + c);
    mesh.AssignPoints(std::move(pts));

    // point attributes
    for (int k = 0; k < num_attrs; ++k) {
        NDArray a(DType::Float64, {static_cast<std::size_t>(npoints)});
        for (std::int64_t i = 0; i < npoints; ++i)
            a.As<double>()[i] = at(i, 4 + k);
        mesh.AddPointData("tetgen:attr" + std::to_string(k + 1), std::move(a));
    }
    // boundary markers: tetgen:ref, tetgen:ref2, ...
    for (int k = 0; k < num_bmarkers; ++k) {
        std::string name = "tetgen:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a(DType::Float64, {static_cast<std::size_t>(npoints)});
        for (std::int64_t i = 0; i < npoints; ++i)
            a.As<double>()[i] = at(i, 4 + num_attrs + k);
        mesh.AddPointData(std::move(name), std::move(a));
    }

    // ---- elements ----
    const detail::FileSource ele_source =
        detail::open_source(ele_path, "Could not open file: " + ele_path);
    const Parsed ef = parse_text(ele_source.View(), ele_path);
    if (ef.mHeader.size() < 3)
        throw ReadError("TetGen: malformed .ele header");
    std::int64_t num_tets = detail::strtoll_token(ef.mHeader[0]);
    int npt = static_cast<int>(detail::strtoll_token(ef.mHeader[1]));
    int ele_attrs = static_cast<int>(detail::strtoll_token(ef.mHeader[2]));
    if (npt != 4)
        throw ReadError("TetGen: only 4-node tetrahedra supported");

    const int ecol = 5 + ele_attrs;
    if (static_cast<std::int64_t>(ef.mData.size()) != num_tets * ecol)
        throw ReadError("TetGen: .ele data size mismatch");

    auto eat = [&](std::int64_t r, int c) -> std::int64_t {
        return detail::strtoll_token(ef.mData[r * ecol + c]);
    };

    NDArray cells(DType::Int64, {static_cast<std::size_t>(num_tets), 4});
    std::int64_t* cp = cells.As<std::int64_t>();
    for (std::int64_t i = 0; i < num_tets; ++i)
        for (int c = 0; c < 4; ++c)
            cp[i * 4 + c] = eat(i, 1 + c) - node_index_base;
    mesh.AddCellBlock("tetra", std::move(cells));

    // region attributes: tetgen:ref, tetgen:ref2, ...
    for (int k = 0; k < ele_attrs; ++k) {
        std::string name = "tetgen:ref" + (k == 0 ? std::string() : std::to_string(k + 1));
        NDArray a(DType::Int64, {static_cast<std::size_t>(num_tets)});
        for (std::int64_t i = 0; i < num_tets; ++i)
            a.As<std::int64_t>()[i] = eat(i, 5 + k);
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(a));
        mesh.AddCellData(std::move(name), std::move(blocks));
    }

    return mesh;
}

namespace {

// Write a marker/ref value: integral values as integers, else %.16e.
void write_value(std::ostream& rOs, double v) {
    double r = std::nearbyint(v);
    if (v == r && std::fabs(v) < 9.2e18) {
        rOs << static_cast<std::int64_t>(r);
    } else {
        char buf[40];
        detail::snprintf_c(buf, sizeof(buf), "%.16e", v);
        rOs << buf;
    }
}

}  // namespace

void write_tetgen(const std::string& rPath, const Mesh& rMesh) {
    bool ok = false;
    auto paths = node_ele_paths(rPath, ok);
    if (!ok)
        throw WriteError("TetGen: must specify a .node or .ele file");
    const std::string& node_path = paths.first;
    const std::string& ele_path = paths.second;

    const NDArray& points = rMesh.Points();
    const std::size_t ncols = rMesh.PointDim();
    if (ncols != 3)
        throw WriteError("TetGen: can only write 3D points");

    const std::int64_t npoints = static_cast<std::int64_t>(rMesh.NumPoints());

    // ---- node file ----
    {
        auto fh = detail::make_classic_ofstream(node_path, std::ios::binary);
        if (!fh)
            throw WriteError("Could not open file for writing: " + node_path);

        // Split point_data into one ref key and the remaining attribute keys,
        // mirroring meshioplusplus.tetgen.write.
        std::vector<std::string> attr_keys =
            rMesh.PointDataNames();  // sorted: deterministic column order
        std::vector<std::string> ref_keys;
        if (!attr_keys.empty()) {
            for (const auto& k : attr_keys)
                if (k.find(":ref") != std::string::npos) {
                    ref_keys.push_back(k);
                    break;
                }
            if (!ref_keys.empty()) {
                attr_keys.erase(std::remove(attr_keys.begin(), attr_keys.end(), ref_keys[0]),
                                attr_keys.end());
            } else {
                ref_keys.push_back(attr_keys.front());
                attr_keys.erase(attr_keys.begin());
            }
        }
        const std::size_t nattr = attr_keys.size();
        const std::size_t nref = ref_keys.size();

        fh << detail::provenance_render_lines(detail::SlotTier::Block, "# ");
        if (nattr + nref > 0) {
            fh << "# attribute and marker names: ";
            bool first = true;
            for (const auto& k : attr_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            for (const auto& k : ref_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            fh << "\n";
        }
        fh << npoints << " 3 " << nattr << " " << nref << "\n";

        char fbuf[40];
        const detail::DoubleView point_values(points);
        std::vector<std::optional<detail::DoubleView>> attrs(attr_keys.size());
        for (std::size_t a = 0; a < attr_keys.size(); ++a)
            attrs[a].emplace(rMesh.PointData(attr_keys[a]));
        std::vector<std::optional<detail::DoubleView>> refs(ref_keys.size());
        for (std::size_t a = 0; a < ref_keys.size(); ++a)
            refs[a].emplace(rMesh.PointData(ref_keys[a]));
        for (std::int64_t i = 0; i < npoints; ++i) {
            fh << i;
            for (int c = 0; c < 3; ++c) {
                detail::snprintf_c(fbuf, sizeof(fbuf), "%.16e", point_values[i * 3 + c]);
                fh << " " << fbuf;
            }
            for (const auto& attr : attrs) {
                detail::snprintf_c(fbuf, sizeof(fbuf), "%.16e", (*attr)[i]);
                fh << " " << fbuf;
            }
            for (const auto& ref : refs) {
                fh << " ";
                write_value(fh, (*ref)[i]);
            }
            fh << "\n";
        }
    }

    // ---- ele file ----
    {
        auto fh = detail::make_classic_ofstream(ele_path, std::ios::binary);
        if (!fh)
            throw WriteError("Could not open file for writing: " + ele_path);

        // Cell-data attribute keys, with the first ":ref" key moved to front.
        std::vector<std::string> attr_keys =
            rMesh.CellDataNames();  // sorted: deterministic column order
        if (!attr_keys.empty()) {
            std::string ref;
            for (const auto& k : attr_keys)
                if (k.find(":ref") != std::string::npos) {
                    ref = k;
                    break;
                }
            if (!ref.empty()) {
                attr_keys.erase(std::remove(attr_keys.begin(), attr_keys.end(), ref),
                                attr_keys.end());
                attr_keys.insert(attr_keys.begin(), ref);
            }
        }
        const std::size_t nattr = attr_keys.size();

        fh << detail::provenance_render_lines(detail::SlotTier::Block, "# ");
        if (nattr > 0) {
            fh << "# attribute names: ";
            bool first = true;
            for (const auto& k : attr_keys) {
                fh << (first ? "" : ", ") << k;
                first = false;
            }
            fh << "\n";
        }

        for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
            const auto cb = rMesh.Cells(ci);
            if (cb.Type() != "tetra")
                continue;
            const NDArray& conn_array = cb.Conn();
            std::int64_t n = detail::rows(conn_array);
            const detail::Int64View conn(conn_array);
            // A key without data for this block writes 0, so its view stays empty.
            std::vector<std::optional<detail::Int64View>> attrs(attr_keys.size());
            for (std::size_t a = 0; a < attr_keys.size(); ++a)
                if (ci < rMesh.CellDataNumBlocks(attr_keys[a]))
                    attrs[a].emplace(rMesh.CellData(attr_keys[a], ci));
            fh << n << " 4 " << nattr << "\n";
            for (std::int64_t i = 0; i < n; ++i) {
                fh << i;
                for (int c = 0; c < 4; ++c)
                    fh << " " << conn[i * 4 + c];
                for (const auto& attr : attrs) {
                    if (attr)
                        fh << " " << (*attr)[i];
                    else
                        fh << " 0";
                }
                fh << "\n";
            }
        }
    }
}

}  // namespace meshioplusplus
