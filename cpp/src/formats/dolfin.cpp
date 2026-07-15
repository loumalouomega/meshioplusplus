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
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

std::pair<std::string, int> dolfin_to_meshio(const std::string& ct) {
    if (ct == "triangle") return {"triangle", 3};
    if (ct == "tetrahedron") return {"tetra", 4};
    throw ReadError("DOLFIN: unsupported cell type '" + ct + "'");
}

const char* meshio_to_dolfin(const std::string& t) {
    if (t == "triangle") return "triangle";
    if (t == "tetra") return "tetrahedron";
    throw WriteError("DOLFIN XML only supports triangles and tetrahedra");
}

}  // namespace

Mesh read_dolfin(const std::string& path) {
    pugi::xml_document doc;
    if (!doc.load_file(path.c_str())) throw ReadError("DOLFIN: could not parse " + path);

    pugi::xml_node dolfin = doc.child("dolfin");
    if (!dolfin) throw ReadError("DOLFIN: missing <dolfin> root");
    pugi::xml_node mesh_node = dolfin.child("mesh");
    if (!mesh_node) throw ReadError("DOLFIN: missing <mesh>");

    int dim = mesh_node.attribute("dim").as_int();
    auto [cell_type, npc] = dolfin_to_meshio(mesh_node.attribute("celltype").value());

    Mesh mesh;

    // Vertices (placed by index).
    pugi::xml_node verts = mesh_node.child("vertices");
    std::size_t nverts = verts.attribute("size").as_uint();
    mesh.points = NDArray(DType::Float64, {nverts, static_cast<std::size_t>(dim)});
    double* pp = mesh.points.as<double>();
    const char* coord[3] = {"x", "y", "z"};
    for (pugi::xml_node v : verts.children("vertex")) {
        std::size_t k = v.attribute("index").as_uint();
        for (int c = 0; c < dim; ++c)
            pp[k * dim + c] = v.attribute(coord[c]).as_double();
    }

    // Cells (single block, placed by index).
    pugi::xml_node cells = mesh_node.child("cells");
    std::size_t ncells = cells.attribute("size").as_uint();
    NDArray data(DType::Int64, {ncells, static_cast<std::size_t>(npc)});
    std::int64_t* dp = data.as<std::int64_t>();
    for (pugi::xml_node c : cells.children()) {
        std::size_t k = c.attribute("index").as_uint();
        for (int j = 0; j < npc; ++j) {
            char tag[8];
            std::snprintf(tag, sizeof(tag), "v%d", j);
            dp[k * npc + j] = c.attribute(tag).as_llong();
        }
    }
    mesh.cells.emplace_back(cell_type, std::move(data));

    // Cell data: sibling files "<stem>_<name>.xml".
    fs::path p(path);
    fs::path dir = p.has_parent_path() ? p.parent_path() : fs::path(".");
    std::string stem = p.stem().string();
    std::string prefix = stem + "_";
    if (fs::exists(dir)) {
        for (const auto& entry : fs::directory_iterator(dir)) {
            std::string fname = entry.path().filename().string();
            if (fname.size() <= prefix.size() + 4) continue;
            if (fname.compare(0, prefix.size(), prefix) != 0) continue;
            if (fname.compare(fname.size() - 4, 4, ".xml") != 0) continue;
            std::string name = fname.substr(prefix.size(), fname.size() - prefix.size() - 4);
            if (name.empty() || name.find('.') != std::string::npos) continue;  // [^.]+

            pugi::xml_document fdoc;
            if (!fdoc.load_file(entry.path().string().c_str())) continue;
            pugi::xml_node mf = fdoc.child("dolfin").child("mesh_function");
            if (!mf) continue;
            std::string type = mf.attribute("type").value();
            std::size_t size = mf.attribute("size").as_uint();
            DType dt = (type == "float") ? DType::Float64 : DType::Int64;
            NDArray arr(dt, {size});
            for (pugi::xml_node e : mf.children("entity")) {
                std::size_t idx = e.attribute("index").as_uint();
                if (dt == DType::Float64)
                    arr.as<double>()[idx] = e.attribute("value").as_double();
                else
                    arr.as<std::int64_t>()[idx] = e.attribute("value").as_llong();
            }
            std::vector<NDArray> blocks;
            blocks.push_back(std::move(arr));
            mesh.cell_data.emplace(name, std::move(blocks));
        }
    }

    return mesh;
}

void write_dolfin(const std::string& path, const Mesh& mesh) {
    // Pick the single supported cell type to write.
    std::string cell_type;
    for (const auto& cb : mesh.cells)
        if (cb.type == "tetra") { cell_type = "tetra"; break; }
    if (cell_type.empty())
        for (const auto& cb : mesh.cells)
            if (cb.type == "triangle") { cell_type = "triangle"; break; }
    if (cell_type.empty())
        throw WriteError("DOLFIN XML only supports triangles and tetrahedra");

    const std::size_t dim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    if (dim != 2 && dim != 3) throw WriteError("DOLFIN: can only write dimension 2 or 3");

    std::ofstream f(path, std::ios::binary);
    if (!f) throw WriteError("Could not open file for writing: " + path);

    f << "<dolfin nsmap=\"{'dolfin': 'https://fenicsproject.org/'}\">\n";
    f << "  <mesh celltype=\"" << meshio_to_dolfin(cell_type) << "\" dim=\"" << dim << "\">\n";

    const std::size_t npts = mesh.num_points();
    f << "    <vertices size=\"" << npts << "\">\n";
    char buf[32];
    const char* coord[3] = {"x", "y", "z"};
    for (std::size_t i = 0; i < npts; ++i) {
        f << "      <vertex index=\"" << i << "\"";
        for (std::size_t c = 0; c < dim; ++c) {
            std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(mesh.points, i * dim + c));
            f << " " << coord[c] << "=\"" << buf << "\"";
        }
        f << " />\n";
    }
    f << "    </vertices>\n";

    std::size_t num_cells = 0;
    for (const auto& cb : mesh.cells)
        if (cb.type == cell_type) num_cells += cb.num_cells();

    f << "    <cells size=\"" << num_cells << "\">\n";
    const char* ts = meshio_to_dolfin(cell_type);
    std::size_t idx = 0;
    for (const auto& cb : mesh.cells) {
        if (cb.type != cell_type) continue;
        std::size_t ncols = detail::cols(cb.data);
        std::size_t n = cb.num_cells();
        for (std::size_t r = 0; r < n; ++r) {
            f << "      <" << ts << " index=\"" << idx << "\"";
            for (std::size_t j = 0; j < ncols; ++j)
                f << " v" << j << "=\"" << detail::read_int(cb.data, r * ncols + j) << "\"";
            f << " />\n";
            ++idx;
        }
    }
    f << "    </cells>\n";
    f << "  </mesh>\n";
    f << "</dolfin>";

    // Cell data -> sibling files "<stem>_<name>.xml".
    bool z_all_zero = true;
    if (dim == 3) {
        for (std::size_t i = 0; i < npts; ++i)
            if (detail::read_double(mesh.points, i * 3 + 2) != 0.0) { z_all_zero = false; break; }
    }
    int data_dim = (dim == 2 || z_all_zero) ? 2 : 3;

    fs::path p(path);
    std::string base = (p.parent_path() / p.stem()).string();
    for (const auto& kv : mesh.cell_data) {
        const std::string& name = kv.first;
        for (const NDArray& arr : kv.second) {
            std::string fn = base + "_" + name + ".xml";
            std::ofstream cf(fn, std::ios::binary);
            if (!cf) throw WriteError("Could not open file for writing: " + fn);
            bool is_float = detail::is_float_dtype(arr.dtype());
            const char* type = is_float ? "float" : "int";
            std::size_t sz = arr.shape().empty() ? 0 : arr.shape()[0];
            cf << "<dolfin><mesh_function type=\"" << type << "\" dim=\"" << data_dim
               << "\" size=\"" << sz << "\">";
            for (std::size_t k = 0; k < sz; ++k) {
                cf << "<entity index=\"" << k << "\" value=\"";
                if (is_float) {
                    std::snprintf(buf, sizeof(buf), "%.17g", detail::read_double(arr, k));
                    cf << buf;
                } else {
                    cf << detail::read_int(arr, k);
                }
                cf << "\" />";
            }
            cf << "</mesh_function></dolfin>";
        }
    }
}

}  // namespace meshioplusplus
