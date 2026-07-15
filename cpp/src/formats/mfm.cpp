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
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// meshio linear type -> (lnv, lne, lnf)
const std::vector<std::pair<std::string, std::array<int, 3>>>& topology() {
    static const std::vector<std::pair<std::string, std::array<int, 3>>> m = {
        {"line", {2, 1, 0}},   {"triangle", {3, 3, 1}}, {"quad", {4, 4, 1}},
        {"tetra", {4, 6, 4}},  {"hexahedron", {8, 12, 6}}, {"wedge", {6, 9, 5}}};
    return m;
}

std::string type_from_dims(int lnv, int lne, int lnf, int lnn) {
    for (const auto& kv : topology())
        if (kv.second[0] == lnv && kv.second[1] == lne && kv.second[2] == lnf &&
            kv.second[0] == lnn)
            return kv.first;
    throw ReadError("MFM: unsupported (non-linear) element");
}

}  // namespace

Mesh read_mfm(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);

    // First non-empty line: header.
    std::string line;
    std::vector<long long> header;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        long long v;
        while (iss >> v) header.push_back(v);
        if (!header.empty()) break;
    }
    if (header.size() < 8) throw ReadError("MFM: expected a header of 8 integers");
    const long long nel = header[0], nnod = header[1], nver = header[2];
    const int dim = static_cast<int>(header[3]);
    const int lnn = static_cast<int>(header[4]), lnv = static_cast<int>(header[5]);
    const int lne = static_cast<int>(header[6]), lnf = static_cast<int>(header[7]);

    std::string cell_type = type_from_dims(lnv, lne, lnf, lnn);
    if (lnn != lnv || nnod != nver)
        throw ReadError("MFM: only linear (P1) elements are supported");

    // Remaining tokens.
    std::vector<std::string> tok((std::istream_iterator<std::string>(in)),
                                 std::istream_iterator<std::string>());
    std::size_t pos = 0;
    auto need = [&](std::size_t n) {
        if (pos + n > tok.size()) throw ReadError("MFM: unexpected end of file");
    };

    NDArray data(DType::Int64, {static_cast<std::size_t>(nel),
                                static_cast<std::size_t>(lnv)});
    need(static_cast<std::size_t>(nel) * lnv);
    for (long long i = 0; i < nel * lnv; ++i)
        data.as<std::int64_t>()[i] = std::strtoll(tok[pos++].c_str(), nullptr, 10) - 1;

    // reference arrays (discarded): nrc (dim==3), nra (dim>=2), nrv
    if (dim == 3) { need(static_cast<std::size_t>(nel) * lnf); pos += nel * lnf; }
    if (dim >= 2) { need(static_cast<std::size_t>(nel) * lne); pos += nel * lne; }
    need(static_cast<std::size_t>(nel) * lnv); pos += nel * lnv;

    Mesh mesh;
    mesh.points = NDArray(DType::Float64,
                          {static_cast<std::size_t>(nver), static_cast<std::size_t>(dim)});
    need(static_cast<std::size_t>(nver) * dim);
    for (long long i = 0; i < nver * dim; ++i)
        mesh.points.as<double>()[i] = std::strtod(tok[pos++].c_str(), nullptr);

    NDArray ref(DType::Int64, {static_cast<std::size_t>(nel)});
    need(static_cast<std::size_t>(nel));
    for (long long i = 0; i < nel; ++i)
        ref.as<std::int64_t>()[i] = std::strtoll(tok[pos++].c_str(), nullptr, 10);

    mesh.cells.emplace_back(cell_type, std::move(data));
    std::vector<NDArray> refs;
    refs.push_back(std::move(ref));
    mesh.cell_data.emplace("mfm:ref", std::move(refs));
    return mesh;
}

void write_mfm(const std::string& path, const Mesh& mesh, const std::string& float_fmt) {
    // Single element type only.
    std::string cell_type;
    for (const auto& cb : mesh.cells) {
        if (cell_type.empty())
            cell_type = cb.type;
        else if (cb.type != cell_type)
            throw WriteError("MFM can only write a single element type");
    }
    if (cell_type.empty()) throw WriteError("MFM: empty mesh");

    const std::array<int, 3>* topo = nullptr;
    for (const auto& kv : topology())
        if (kv.first == cell_type) { topo = &kv.second; break; }
    if (!topo) throw WriteError("MFM does not support '" + cell_type + "' cells");
    const int lnv = (*topo)[0], lne = (*topo)[1], lnf = (*topo)[2], lnn = lnv;

    std::size_t nel = 0;
    for (const auto& cb : mesh.cells) nel += cb.num_cells();
    const std::size_t nver = mesh.num_points();
    const int dim = mesh.points.shape().size() >= 2
                        ? static_cast<int>(mesh.points.shape()[1])
                        : 0;

    // subdomain
    std::vector<std::int64_t> nsd(nel, 1);
    auto it = mesh.cell_data.find("mfm:ref");
    if (it != mesh.cell_data.end()) {
        std::size_t p = 0;
        for (const auto& blk : it->second)
            for (std::size_t i = 0; i < blk.size() && p < nel; ++i)
                nsd[p++] = detail::read_int(blk, i);
    }

    std::ofstream f(path, std::ios::binary);
    if (!f) throw WriteError("Could not open file for writing: " + path);
    f << nel << " " << nver << " " << nver << " " << dim << " " << lnn << " " << lnv
      << " " << lne << " " << lnf << "\n";

    // connectivity (1-based)
    for (const auto& cb : mesh.cells) {
        std::size_t n = cb.num_cells();
        std::size_t k = detail::cols(cb.data);
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t j = 0; j < k; ++j)
                f << (detail::read_int(cb.data, r * k + j) + 1) << (j + 1 == k ? '\n' : ' ');
        }
    }
    // zero reference arrays
    auto zeros = [&](int cols) {
        for (std::size_t r = 0; r < nel; ++r)
            for (int j = 0; j < cols; ++j) f << 0 << (j + 1 == cols ? '\n' : ' ');
    };
    if (dim == 3) zeros(lnf);
    if (dim >= 2) zeros(lne);
    zeros(lnv);
    // coordinates
    std::string fmt = "%" + float_fmt;
    char buf[64];
    for (std::size_t i = 0; i < nver; ++i)
        for (int c = 0; c < dim; ++c) {
            std::snprintf(buf, sizeof(buf), fmt.c_str(),
                          detail::read_double(mesh.points, i * dim + c));
            f << buf << (c + 1 == dim ? '\n' : ' ');
        }
    // subdomain
    for (std::size_t i = 0; i < nel; ++i) f << nsd[i] << "\n";
}

}  // namespace meshioplusplus
