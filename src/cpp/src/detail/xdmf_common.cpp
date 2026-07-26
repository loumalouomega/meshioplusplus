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
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <unordered_map>

// Project includes
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/detail/hdf5_util.hpp"
#endif

namespace meshioplusplus {
namespace xdmfcommon {

const char* meshio_to_xdmf(const std::string& rT) {
    static const std::unordered_map<std::string, const char*> m = {
        {"vertex", "Polyvertex"},
        {"line", "Polyline"},
        {"line3", "Edge_3"},
        {"quad", "Quadrilateral"},
        {"quad8", "Quadrilateral_8"},
        {"quad9", "Quadrilateral_9"},
        {"pyramid", "Pyramid"},
        {"pyramid13", "Pyramid_13"},
        {"tetra", "Tetrahedron"},
        {"triangle", "Triangle"},
        {"triangle6", "Triangle_6"},
        {"tetra10", "Tetrahedron_10"},
        {"wedge", "Wedge"},
        {"wedge15", "Wedge_15"},
        {"wedge18", "Wedge_18"},
        {"hexahedron", "Hexahedron"},
        {"hexahedron20", "Hexahedron_20"},
        {"hexahedron24", "Hexahedron_24"},
        {"hexahedron27", "Hexahedron_27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw WriteError("XDMF: unsupported cell type " + rT);
    return it->second;
}

std::string xdmf_to_meshio(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"Polyvertex", "vertex"},
        {"Polyline", "line"},
        {"Edge_3", "line3"},
        {"Quadrilateral", "quad"},
        {"Quadrilateral_8", "quad8"},
        {"Quad_8", "quad8"},
        {"Quadrilateral_9", "quad9"},
        {"Quad_9", "quad9"},
        {"Pyramid", "pyramid"},
        {"Pyramid_13", "pyramid13"},
        {"Tetrahedron", "tetra"},
        {"Triangle", "triangle"},
        {"Triangle_6", "triangle6"},
        {"Tri_6", "triangle6"},
        {"Tetrahedron_10", "tetra10"},
        {"Tet_10", "tetra10"},
        {"Wedge", "wedge"},
        {"Wedge_15", "wedge15"},
        {"Wedge_18", "wedge18"},
        {"Hexahedron", "hexahedron"},
        {"Hexahedron_20", "hexahedron20"},
        {"Hex_20", "hexahedron20"},
        {"Hexahedron_24", "hexahedron24"},
        {"Hex_24", "hexahedron24"},
        {"Hexahedron_27", "hexahedron27"},
        {"Hex_27", "hexahedron27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw ReadError("XDMF: unsupported topology type " + rT);
    return it->second;
}

NDArray concat_cell_data(const Mesh& rMesh, const std::string& rName) {
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::size_t total_rows = 0;
    std::vector<std::size_t> shape = rMesh.CellData(rName, 0).Shape();
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto& bshape = rMesh.CellData(rName, b).Shape();
        total_rows += bshape.empty() ? 0 : bshape[0];
    }
    shape[0] = total_rows;
    NDArray out(rMesh.CellData(rName, 0).Dtype(), shape);
    std::size_t off = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& blk = rMesh.CellData(rName, b);
        std::memcpy(out.Data() + off, blk.Data(), blk.Nbytes());
        off += blk.Nbytes();
    }
    return out;
}

std::vector<NDArray> split_raw_cell_data(const NDArray& rRaw,
                                         const std::vector<std::size_t>& rSizes) {
    std::size_t ncols = rRaw.Ndim() >= 2 ? rRaw.Shape()[1] : 1;
    std::size_t off = 0;
    std::vector<NDArray> blocks;
    for (std::size_t bs : rSizes) {
        std::vector<std::size_t> bshape = rRaw.Shape();
        if (!bshape.empty())
            bshape[0] = bs;
        NDArray b(rRaw.Dtype(), bshape);
        std::size_t elems = bs * ncols;
        std::memcpy(b.Data(), rRaw.Data() + off * ncols * dtype_size(rRaw.Dtype()),
                    elems * dtype_size(rRaw.Dtype()));
        off += bs;
        blocks.push_back(std::move(b));
    }
    return blocks;
}

// ---- write-path helpers ----

std::pair<const char*, const char*> numpy_to_xdmf_dtype(DType dt) {
    switch (dt) {
        case DType::Int8:
            return {"Int", "1"};
        case DType::Int16:
            return {"Int", "2"};
        case DType::Int32:
            return {"Int", "4"};
        case DType::Int64:
            return {"Int", "8"};
        case DType::UInt8:
            return {"UInt", "1"};
        case DType::UInt16:
            return {"UInt", "2"};
        case DType::UInt32:
            return {"UInt", "4"};
        case DType::UInt64:
            return {"UInt", "8"};
        case DType::Float32:
            return {"Float", "4"};
        case DType::Float64:
            return {"Float", "8"};
    }
    return {"Float", "8"};
}

std::string attribute_type(const std::vector<std::size_t>& rShape) {
    if (rShape.size() == 1 || (rShape.size() == 2 && rShape[1] == 1))
        return "Scalar";
    if (rShape.size() == 2 && (rShape[1] == 2 || rShape[1] == 3))
        return "Vector";
    if ((rShape.size() == 2 && rShape[1] == 9) ||
        (rShape.size() == 3 && rShape[1] == 3 && rShape[2] == 3))
        return "Tensor";
    if (rShape.size() == 2 && rShape[1] == 6)
        return "Tensor6";
    return "Matrix";
}

int meshio_to_xdmf_index(const std::string& rT) {
    static const std::unordered_map<std::string, int> m = {
        {"vertex", 0x1},        {"line", 0x2},          {"triangle", 0x4},     {"quad", 0x5},
        {"tetra", 0x6},         {"pyramid", 0x7},       {"wedge", 0x8},        {"hexahedron", 0x9},
        {"line3", 0x22},        {"quad9", 0x23},        {"triangle6", 0x24},   {"quad8", 0x25},
        {"tetra10", 0x26},      {"pyramid13", 0x27},    {"wedge15", 0x28},     {"wedge18", 0x29},
        {"hexahedron20", 0x30}, {"hexahedron24", 0x31}, {"hexahedron27", 0x32}};
    auto it = m.find(rT);
    if (it == m.end())
        throw WriteError("XDMF: cannot mix cell type " + rT);
    return it->second;
}

std::string dims_string(const NDArray& rArr) {
    std::string dims;
    for (std::size_t i = 0; i < rArr.Shape().size(); ++i) {
        if (i)
            dims += " ";
        dims += std::to_string(rArr.Shape()[i]);
    }
    return dims;
}

NDArray pack_mixed_topology(const Mesh& rMesh, std::size_t& rTotalCells) {
    std::size_t total_cells = 0, total_len = 0;
    for (const auto cb : rMesh.CellRange()) {
        std::size_t nc = cb.NumCells();
        std::size_t npc = detail::cols(cb.Conn());
        std::size_t prefix = (cb.Type() == "vertex" || cb.Type() == "line") ? 2 : 1;
        total_cells += nc;
        total_len += nc * (prefix + npc);
    }
    NDArray cd(DType::Int64, {total_len});
    std::int64_t* cp = cd.As<std::int64_t>();
    std::size_t pos = 0;
    for (const auto cb : rMesh.CellRange()) {
        std::size_t nc = cb.NumCells();
        const NDArray& conn = cb.Conn();
        std::size_t npc = detail::cols(conn);
        int idx = meshio_to_xdmf_index(cb.Type());
        std::size_t prefix = (cb.Type() == "vertex" || cb.Type() == "line") ? 2 : 1;
        for (std::size_t r = 0; r < nc; ++r) {
            for (std::size_t pq = 0; pq < prefix; ++pq)
                cp[pos++] = idx;
            for (std::size_t j = 0; j < npc; ++j)
                cp[pos++] = detail::read_int(conn, r * npc + j);
        }
    }
    rTotalCells = total_cells;
    return cd;
}

struct DataItemStore::Impl {
    std::string mDataFormat;
    std::string mBase;
    int mCounter = 0;
    int mGzipLevel = -1;
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    h5::Hid mH5File;
    std::string mH5Basename;
#endif
};

DataItemStore::DataItemStore(const std::string& rDataFormat, const std::string& rBase,
                             int GzipLevel)
    : mImpl(std::make_unique<Impl>()) {
    mImpl->mDataFormat = rDataFormat;
    mImpl->mBase = rBase;
    mImpl->mGzipLevel = GzipLevel;
}

DataItemStore::~DataItemStore() = default;

const std::string& DataItemStore::DataFormat() const {
    return mImpl->mDataFormat;
}

void DataItemStore::Close() {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    mImpl->mH5File.Reset();
#endif
}

void DataItemStore::Flush() {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    h5::flush_file(mImpl->mH5File);
#endif
    // "Binary" opens, writes and closes one file per array, and "XML" keeps its
    // numbers in the caller's document -- neither has anything buffered here.
}

int DataItemStore::Counter() const {
    return mImpl->mCounter;
}

void DataItemStore::SetCounter(int Value) {
    mImpl->mCounter = Value;
}

void DataItemStore::OpenExisting() {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    if (mImpl->mDataFormat != "HDF")
        return;
    const std::string h5_path = mImpl->mBase + ".h5";
    std::error_code ec;
    if (!std::filesystem::exists(h5_path, ec))
        return;  // nothing to continue; Store() will create it as usual
    mImpl->mH5File = h5::open_file_rw(h5_path);
    const std::size_t slash = h5_path.find_last_of("/\\");
    mImpl->mH5Basename = slash == std::string::npos ? h5_path : h5_path.substr(slash + 1);
    // Resume past every dataN already there. Scanning the container rather than
    // trusting the document is deliberate: a mis-resumed counter would silently
    // overwrite data0 rather than fail, and the file is the authority on what
    // it holds.
    int next = 0;
    for (const std::string& r_name : h5::link_names(mImpl->mH5File)) {
        if (r_name.rfind("data", 0) != 0)
            continue;
        const std::string digits = r_name.substr(4);
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        next = std::max(next, std::stoi(digits) + 1);
    }
    mImpl->mCounter = std::max(mImpl->mCounter, next);
#endif
}

std::string DataItemStore::Store(const NDArray& rArr) {
    const std::size_t rows = rArr.Shape().empty() ? 0 : rArr.Shape()[0];
    const std::size_t cols = rows ? rArr.Size() / rows : 0;

    if (mImpl->mDataFormat == "Binary") {
        std::string fn = mImpl->mBase + std::to_string(mImpl->mCounter++) + ".bin";
        std::ofstream bf(fn, std::ios::binary);
        if (!bf)
            throw WriteError("XDMF: could not write " + fn);
        bf.write(reinterpret_cast<const char*>(rArr.Data()),
                 static_cast<std::streamsize>(rArr.Nbytes()));
        return fn;
    }
    if (mImpl->mDataFormat == "HDF") {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        if (!mImpl->mH5File.Valid()) {
            std::string h5_path = mImpl->mBase + ".h5";
            mImpl->mH5File = h5::create_file(h5_path);
            std::size_t slash = h5_path.find_last_of("/\\");
            mImpl->mH5Basename = slash == std::string::npos ? h5_path : h5_path.substr(slash + 1);
        }
        std::string name = "data" + std::to_string(mImpl->mCounter++);
        h5::write_dataset(mImpl->mH5File, name, rArr, mImpl->mGzipLevel);
        return mImpl->mH5Basename + ":/" + name;
#else
        throw WriteError(
            "XDMF: HDF data format requires an HDF5-enabled build "
            "(-DMESHIOPLUSPLUS_WITH_HDF5=ON)");
#endif
    }
    // XML inline
    std::string text = "\n";
    char buf[40];
    const bool is_float = detail::is_float_dtype(rArr.Dtype());
    const bool f32 = rArr.Dtype() == DType::Float32;
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t cc = 0; cc < cols; ++cc) {
            std::size_t i = r * cols + cc;
            if (is_float) {
                std::snprintf(buf, sizeof(buf), f32 ? "%.7e" : "%.16e",
                              detail::read_double(rArr, i));
            } else {
                std::snprintf(buf, sizeof(buf), "%lld",
                              static_cast<long long>(detail::read_int(rArr, i)));
            }
            if (cc)
                text += " ";
            text += buf;
        }
        text += "\n";
    }
    return text;
}

}  // namespace xdmfcommon
}  // namespace meshioplusplus
