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
// GIFTI surfaces. See formats/gifti.hpp.

// System includes
#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/gifti.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"
#include "pugixml.hpp"

namespace meshioplusplus {

namespace {

DType gii_dtype(const std::string& rName) {
    static const std::map<std::string, DType> table = {
        {"NIFTI_TYPE_UINT8", DType::UInt8},     {"NIFTI_TYPE_INT8", DType::Int8},
        {"NIFTI_TYPE_INT16", DType::Int16},     {"NIFTI_TYPE_UINT16", DType::UInt16},
        {"NIFTI_TYPE_INT32", DType::Int32},     {"NIFTI_TYPE_UINT32", DType::UInt32},
        {"NIFTI_TYPE_INT64", DType::Int64},     {"NIFTI_TYPE_UINT64", DType::UInt64},
        {"NIFTI_TYPE_FLOAT32", DType::Float32}, {"NIFTI_TYPE_FLOAT64", DType::Float64}};
    const auto it = table.find(rName);
    if (it == table.end())
        throw ReadError("meshio++: gifti: DataType '" + rName + "' is not supported");
    return it->second;
}

/// The `<MD>` value named @p rKey in @p Node's `<MetaData>`, or "".
std::string gii_meta(const pugi::xml_node& Node, const std::string& rKey) {
    for (const auto md : Node.child("MetaData").children("MD"))
        if (std::string(md.child("Name").text().get()) == rKey)
            return md.child("Value").text().get();
    return {};
}

std::string gii_trim(std::string Text) {
    const auto b = Text.find_first_not_of(" \t\r\n");
    const auto e = Text.find_last_not_of(" \t\r\n");
    return b == std::string::npos ? std::string() : Text.substr(b, e - b + 1);
}

/// One `DataArray`'s values, row-major, of its own dtype.
NDArray gii_decode(const pugi::xml_node& Array) {
    const DType dt = gii_dtype(Array.attribute("DataType").value());
    const int rank = Array.attribute("Dimensionality").as_int(1);
    if (rank < 1 || rank > 4)
        throw ReadError("meshio++: gifti: Dimensionality must be 1 to 4");
    std::vector<std::size_t> shape;
    std::size_t count = 1;
    for (int d = 0; d < rank; ++d) {
        const long long n = Array.attribute(("Dim" + std::to_string(d)).c_str()).as_llong(-1);
        if (n < 0 || n > (1LL << 40))
            throw ReadError("meshio++: gifti: Dim" + std::to_string(d) + " is missing or invalid");
        shape.push_back(static_cast<std::size_t>(n));
        count *= static_cast<std::size_t>(n);
    }
    const std::string encoding = Array.attribute("Encoding").value();
    const std::string text = Array.child("Data").text().get();
    NDArray out(dt, shape);
    const std::size_t itemsize = out.Nbytes() / std::max<std::size_t>(count, 1);
    if (encoding == "ASCII") {
        const char* p = text.c_str();
        for (std::size_t k = 0; k < count; ++k) {
            while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
                ++p;
            const char* end = nullptr;
            const double v = detail::parse_double(p, end);
            if (end == p)
                throw ReadError(
                    "meshio++: gifti: an ASCII DataArray holds fewer values than its "
                    "dimensions");
            p = end;
            detail::write_double(out, k, v);
        }
    } else if (encoding == "Base64Binary" || encoding == "GZipBase64Binary") {
        std::vector<unsigned char> raw = detail::b64decode(text.data(), text.size());
        std::string bytes(raw.begin(), raw.end());
        if (encoding == "GZipBase64Binary")
            bytes = detail::zlib_inflate(bytes, 15, nullptr, "gifti");
        if (bytes.size() != out.Nbytes())
            throw ReadError("meshio++: gifti: a DataArray holds " + std::to_string(bytes.size()) +
                            " bytes, not the " + std::to_string(out.Nbytes()) +
                            " its dimensions need");
        if (!bytes.empty())
            std::memcpy(out.Data(), bytes.data(), bytes.size());
        const bool big = std::string(Array.attribute("Endian").value()) == "BigEndian";
        if (big != (std::endian::native == std::endian::big) && itemsize > 1) {
            auto* b = reinterpret_cast<unsigned char*>(out.Data());
            for (std::size_t k = 0; k < count; ++k)
                std::reverse(b + k * itemsize, b + (k + 1) * itemsize);
        }
    } else {
        throw ReadError("meshio++: gifti: Encoding '" + encoding +
                        "' is not read (ExternalFileBinary needs its sibling file)");
    }
    if (std::string(Array.attribute("ArrayIndexingOrder").value()) == "ColumnMajorOrder" &&
        rank == 2 && shape[0] > 0 && shape[1] > 1) {
        // Column-major storage: element (i, j) sits at j * n0 + i.
        NDArray t(dt, shape);
        const auto* src = reinterpret_cast<const unsigned char*>(out.Data());
        auto* dst = reinterpret_cast<unsigned char*>(t.Data());
        for (std::size_t i = 0; i < shape[0]; ++i)
            for (std::size_t j = 0; j < shape[1]; ++j)
                std::memcpy(dst + (i * shape[1] + j) * itemsize,
                            src + (j * shape[0] + i) * itemsize, itemsize);
        out = std::move(t);
    }
    return out;
}

NDArray gii_text(const std::string& rText) {
    NDArray out(DType::UInt8, {rText.size()});
    if (!rText.empty())
        std::memcpy(out.Data(), rText.data(), rText.size());
    return out;
}

}  // namespace

Mesh read_gifti(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: gifti: could not open file: " + rPath);
    pugi::xml_document doc;
    const std::string_view view = source.View();
    if (!doc.load_buffer(view.data(), view.size()))
        throw ReadError("meshio++: gifti: could not parse " + rPath + " as XML");
    const pugi::xml_node root = doc.child("GIFTI");
    if (!root)
        throw ReadError("meshio++: gifti: '" + rPath + "' has no <GIFTI> root");
    Mesh mesh;
    bool have_points = false;
    std::size_t np = 0;
    std::vector<std::pair<std::string, NDArray>> later;
    std::vector<std::pair<std::string, bool>> later_label;
    std::size_t index = 0;
    for (const auto array : root.children("DataArray")) {
        const std::string intent = array.attribute("Intent").value();
        NDArray data = gii_decode(array);
        if (intent == "NIFTI_INTENT_POINTSET" && !have_points) {
            if (data.Shape().size() != 2 || data.Shape()[1] != 3)
                throw ReadError("meshio++: gifti: the POINTSET is not (n, 3)");
            np = data.Shape()[0];
            NDArray pts(DType::Float64, {np, 3});
            for (std::size_t k = 0; k < np * 3; ++k)
                pts.As<double>()[k] = detail::read_double(data, k);
            mesh.AssignPoints(std::move(pts));
            have_points = true;
            const auto xform = array.child("CoordinateSystemTransformMatrix");
            if (xform) {
                NDArray m(DType::Float64, {4, 4});
                const std::string text = xform.child("MatrixData").text().get();
                const char* p = text.c_str();
                for (std::size_t k = 0; k < 16; ++k) {
                    while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')
                        ++p;
                    const char* end = nullptr;
                    m.As<double>()[k] = detail::parse_double(p, end);
                    if (end == p)
                        throw ReadError("meshio++: gifti: MatrixData needs 16 numbers");
                    p = end;
                }
                mesh.AddFieldData("gifti:transform", std::move(m));
                mesh.AddFieldData("gifti:data_space",
                                  gii_text(gii_trim(xform.child("DataSpace").text().get())));
                mesh.AddFieldData("gifti:transformed_space",
                                  gii_text(gii_trim(xform.child("TransformedSpace").text().get())));
            }
        } else if (intent == "NIFTI_INTENT_TRIANGLE") {
            if (data.Shape().size() != 2 || data.Shape()[1] != 3)
                throw ReadError("meshio++: gifti: the TRIANGLE array is not (n, 3)");
            NDArray conn(DType::Int64, data.Shape());
            for (std::size_t k = 0; k < data.Size(); ++k)
                conn.As<std::int64_t>()[k] = detail::read_int(data, k);
            later.emplace_back("\x01triangle", std::move(conn));
        } else {
            std::string name = gii_meta(array, "Name");
            if (name.empty())
                name = intent == "NIFTI_INTENT_LABEL"
                           ? std::string("gifti:label")
                           : "gifti:" + intent + "_" + std::to_string(index);
            later_label.emplace_back(name, intent == "NIFTI_INTENT_LABEL");
            later.emplace_back(name, std::move(data));
        }
        ++index;
    }
    if (!have_points)
        throw ReadError("meshio++: gifti: '" + rPath + "' has no NIFTI_INTENT_POINTSET array");
    std::size_t data_k = 0;
    std::string label_name;
    for (auto& [name, data] : later) {
        if (name == "\x01triangle") {
            for (std::size_t k = 0; k < data.Size(); ++k)
                if (data.As<std::int64_t>()[k] < 0 ||
                    static_cast<std::size_t>(data.As<std::int64_t>()[k]) >= np)
                    throw ReadError("meshio++: gifti: a triangle names a point out of range");
            mesh.AddCellBlock("triangle", std::move(data));
            continue;
        }
        const bool is_label = later_label[data_k++].second;
        if (data.Shape().empty() || data.Shape()[0] != np) {
            log::warn("gifti: array '{}' does not hold one row per point; skipped", name);
            continue;
        }
        if (data.Shape().size() > 2) {
            std::size_t width = data.Size() / np;
            data.Reshape({np, width});
        }
        if (is_label && label_name.empty())
            label_name = name;
        mesh.AddPointData(name, std::move(data));
    }
    // The LabelTable's keys become point regions over the label array.
    const pugi::xml_node table = root.child("LabelTable");
    if (table && !label_name.empty()) {
        const NDArray& labels = mesh.PointData(label_name);
        for (const auto label : table.children("Label")) {
            const long long key =
                label.attribute("Key").as_llong(label.attribute("Index").as_llong());
            std::vector<std::int64_t> ids;
            for (std::size_t k = 0; k < np; ++k)
                if (detail::read_int(labels, k) == key)
                    ids.push_back(static_cast<std::int64_t>(k));
            NDArray entries(DType::Int64, {ids.size()});
            std::copy(ids.begin(), ids.end(), entries.As<std::int64_t>());
            mesh.AddRegion(Region(gii_trim(label.text().get()), RegionKind::Point, -1, key,
                                  std::move(entries)));
        }
    }
    return mesh;
}

namespace {

std::string gii_array(const std::string& rIntent, const std::string& rType,
                      const std::vector<std::size_t>& rShape, const std::string& rBytes,
                      const std::string& rName, const std::string& rExtra) {
    std::string s = "<DataArray Intent=\"" + rIntent + "\" DataType=\"" + rType +
                    "\" ArrayIndexingOrder=\"RowMajorOrder\" Dimensionality=\"" +
                    std::to_string(rShape.size()) + "\"";
    for (std::size_t d = 0; d < rShape.size(); ++d)
        s += " Dim" + std::to_string(d) + "=\"" + std::to_string(rShape[d]) + "\"";
    s += " Encoding=\"Base64Binary\" Endian=\"LittleEndian\" ExternalFileName=\"\" "
         "ExternalFileOffset=\"\">\n";
    s += rName.empty() ? "<MetaData/>\n"
                       : "<MetaData>\n<MD>\n<Name><![CDATA[Name]]></Name>\n<Value><![CDATA[" +
                             rName + "]]></Value>\n</MD>\n</MetaData>\n";
    s += rExtra;
    s += "<Data>" +
         detail::b64encode(reinterpret_cast<const unsigned char*>(rBytes.data()), rBytes.size()) +
         "</Data>\n</DataArray>\n";
    return s;
}

template <class T>
void gii_put(std::string& rOut, T Value) {
    if constexpr (std::endian::native == std::endian::big) {
        unsigned char b[sizeof(T)];
        std::memcpy(b, &Value, sizeof(T));
        std::reverse(b, b + sizeof(T));
        rOut.append(reinterpret_cast<const char*>(b), sizeof(T));
    } else {
        rOut.append(reinterpret_cast<const char*>(&Value), sizeof(T));
    }
}

}  // namespace

void write_gifti(const std::string& rPath, const Mesh& rMesh) {
    if (rMesh.NumCellBlocks() > 1 ||
        (rMesh.NumCellBlocks() == 1 && std::string(rMesh.Cells(0).Type()) != "triangle"))
        throw WriteError("meshio++: gifti: a GIFTI surface is one triangle block");
    const std::size_t np = rMesh.NumPoints(), pd = rMesh.PointDim();
    std::vector<std::string> arrays;
    std::string bytes;
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            gii_put(bytes, static_cast<float>(
                               c < pd ? detail::read_double(rMesh.Points(), i * pd + c) : 0.0));
    std::string xform = "<CoordinateSystemTransformMatrix>\n<DataSpace><![CDATA[";
    const auto text_of = [&](const char* pKey, const char* pDefault) {
        if (!rMesh.HasFieldData(pKey))
            return std::string(pDefault);
        const NDArray& t = rMesh.FieldData(pKey);
        return std::string(reinterpret_cast<const char*>(t.Data()), t.Size());
    };
    xform += text_of("gifti:data_space", "NIFTI_XFORM_UNKNOWN") +
             "]]></DataSpace>\n<TransformedSpace><![CDATA[" +
             text_of("gifti:transformed_space", "NIFTI_XFORM_UNKNOWN") +
             "]]></TransformedSpace>\n<MatrixData>";
    char buf[64];
    for (std::size_t k = 0; k < 16; ++k) {
        double v = (k % 5 == 0) ? 1.0 : 0.0;
        if (rMesh.HasFieldData("gifti:transform") &&
            rMesh.FieldData("gifti:transform").Size() == 16)
            v = detail::read_double(rMesh.FieldData("gifti:transform"), k);
        detail::snprintf_c(buf, sizeof buf, "%.17g", v);
        xform += (k ? " " : "") + std::string(buf);
    }
    xform += "</MatrixData>\n</CoordinateSystemTransformMatrix>\n";
    arrays.push_back(
        gii_array("NIFTI_INTENT_POINTSET", "NIFTI_TYPE_FLOAT32", {np, 3}, bytes, "", xform));
    if (rMesh.NumCellBlocks() == 1) {
        const auto cb = rMesh.Cells(0);
        bytes.clear();
        for (std::size_t k = 0; k < cb.NumCells() * 3; ++k)
            gii_put(bytes, static_cast<std::int32_t>(detail::read_int(cb.Conn(), k)));
        arrays.push_back(gii_array("NIFTI_INTENT_TRIANGLE", "NIFTI_TYPE_INT32", {cb.NumCells(), 3},
                                   bytes, "", ""));
    }
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (a.Shape().empty() || a.Shape()[0] != np)
            continue;
        const bool floating = a.Dtype() == DType::Float32 || a.Dtype() == DType::Float64;
        const bool integer = !floating;
        if (!floating && !integer)
            continue;
        const std::size_t width = np ? a.Size() / np : 0;
        bytes.clear();
        for (std::size_t k = 0; k < a.Size(); ++k) {
            if (floating)
                gii_put(bytes, static_cast<float>(detail::read_double(a, k)));
            else
                gii_put(bytes, static_cast<std::int32_t>(detail::read_int(a, k)));
        }
        const std::string intent = name == "gifti:label" ? "NIFTI_INTENT_LABEL"
                                   : width > 1           ? "NIFTI_INTENT_NONE"
                                   : floating            ? "NIFTI_INTENT_SHAPE"
                                                         : "NIFTI_INTENT_NONE";
        arrays.push_back(gii_array(
            intent, floating ? "NIFTI_TYPE_FLOAT32" : "NIFTI_TYPE_INT32",
            width > 1 ? std::vector<std::size_t>{np, width} : std::vector<std::size_t>{np}, bytes,
            name, ""));
    }
    std::string out =
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<!DOCTYPE GIFTI SYSTEM \"http://www.nitrc.org/frs/download.php/115/"
        "gifti.dtd\">\n<GIFTI Version=\"1.0\" NumberOfDataArrays=\"" +
        std::to_string(arrays.size()) + "\">\n<MetaData/>\n";
    std::string labels;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Point && r.mTag >= 0)
            labels += "<Label Key=\"" + std::to_string(r.mTag) +
                      "\" Red=\"0\" Green=\"0\" Blue=\"0\" Alpha=\"0\"><![CDATA[" + r.mName +
                      "]]></Label>\n";
    }
    out += labels.empty() ? "<LabelTable/>\n" : "<LabelTable>\n" + labels + "</LabelTable>\n";
    for (const std::string& a : arrays)
        out += a;
    out += "</GIFTI>\n";
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: gifti: could not open file for writing: " + rPath);
    os << out;
    if (!os)
        throw WriteError("meshio++: gifti: failed writing " + rPath);
}

}  // namespace meshioplusplus
