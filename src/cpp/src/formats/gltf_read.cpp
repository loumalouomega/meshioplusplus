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
// glTF 2.0 reader. See formats/gltf.hpp for the contract. The Python twin
// (src/python/meshioplusplus/gltf/_gltf_read.py) follows the same steps.

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <ios>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/read_options.hpp"
#include "../detail/json_value.hpp"

namespace meshioplusplus {
namespace {

using detail::JsonValue;
using GltfrMat4 = std::array<double, 16>;  // column-major, as glTF writes `matrix`

constexpr std::uint64_t kGltfrMaxCount = 1ull << 27;
constexpr std::uint64_t kGltfrMaxZeros = 1ull << 26;
constexpr std::uint64_t kGltfrMaxInt = 1ull << 40;
constexpr double kGltfrHalfSqrt2 = 0.70710678118654757;

[[noreturn]] void gltfr_fail(const std::string& rMessage) {
    throw ReadError("gltf: " + rMessage);
}

// ---- the JSON accessors -------------------------------------------------- //

const JsonValue* gltfr_get(const JsonValue& rObject, const char* pKey) {
    return rObject.IsObject() ? rObject.Find(pKey) : nullptr;
}

const JsonValue& gltfr_object(const JsonValue* pValue, const std::string& rWhat) {
    if (pValue == nullptr || !pValue->IsObject())
        gltfr_fail(rWhat + " is not an object");
    return *pValue;
}

std::uint64_t gltfr_uint(const JsonValue& rValue, const std::string& rWhat) {
    if (!rValue.IsNumber() || !rValue.mIsInteger || rValue.mInteger < 0 ||
        static_cast<std::uint64_t>(rValue.mInteger) > kGltfrMaxInt)
        gltfr_fail(rWhat + " is not a non-negative integer");
    return static_cast<std::uint64_t>(rValue.mInteger);
}

std::uint64_t gltfr_uint_or(const JsonValue& rObject, const char* pKey, std::uint64_t Default,
                            const std::string& rWhat) {
    const JsonValue* v = gltfr_get(rObject, pKey);
    return v == nullptr ? Default : gltfr_uint(*v, rWhat + "." + pKey);
}

// The array under @p pKey, or an empty list.
const std::vector<JsonValue>& gltfr_list(const JsonValue& rObject, const char* pKey) {
    static const std::vector<JsonValue> kEmpty;
    const JsonValue* v = gltfr_get(rObject, pKey);
    if (v == nullptr)
        return kEmpty;
    if (!v->IsArray())
        gltfr_fail(std::string("'") + pKey + "' is not an array");
    return v->mItems;
}

std::vector<double> gltfr_numbers(const JsonValue& rValue, std::size_t Count,
                                  const std::string& rWhat) {
    if (!rValue.IsArray() || rValue.mItems.size() != Count)
        gltfr_fail(rWhat + " is not an array of " + std::to_string(Count) + " numbers");
    std::vector<double> out;
    out.reserve(Count);
    for (const JsonValue& v : rValue.mItems) {
        if (!v.IsNumber() || !std::isfinite(v.mNumber))
            gltfr_fail(rWhat + " holds a value that is not a finite number");
        out.push_back(v.mNumber);
    }
    return out;
}

// ---- the container ------------------------------------------------------- //

bool gltfr_has_suffix_ci(const std::string& rText, const char* pSuffix) {
    const std::string s = pSuffix;
    if (rText.size() < s.size())
        return false;
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = rText[rText.size() - s.size() + i];
        if (c >= 'A' && c <= 'Z')
            c = static_cast<char>(c - 'A' + 'a');
        if (c != s[i])
            return false;
    }
    return true;
}

std::string gltfr_read_file(const std::string& rPath) {
    auto f = detail::make_classic_ifstream(rPath, std::ios::in | std::ios::binary);
    if (!f)
        gltfr_fail("cannot open '" + rPath + "'");
    const std::size_t size = detail::file_bytes(rPath);
    std::string data(size, '\0');
    if (size > 0)
        f.read(data.data(), static_cast<std::streamsize>(size));
    data.resize(static_cast<std::size_t>(f.gcount()));
    return data;
}

std::uint32_t gltfr_u32(const std::string& rData, std::size_t Offset) {
    std::uint32_t v = 0;
    for (int k = 3; k >= 0; --k)
        v = (v << 8) | static_cast<unsigned char>(rData[Offset + static_cast<std::size_t>(k)]);
    return v;
}

struct GltfrFile {
    JsonValue mRoot;
    std::string mBin;
    bool mHasBin = false;
    std::string mDirectory;
};

GltfrFile gltfr_load(const std::string& rPath) {
    GltfrFile file;
    const std::size_t cut = rPath.find_last_of("/\\");
    file.mDirectory = cut == std::string::npos ? std::string() : rPath.substr(0, cut + 1);
    const std::string data = gltfr_read_file(rPath);

    std::string json_text;
    if (data.size() >= 4 && std::memcmp(data.data(), "glTF", 4) == 0) {
        if (data.size() < 20)
            gltfr_fail("the GLB is truncated");
        if (gltfr_u32(data, 4) != 2u)
            gltfr_fail("only GLB version 2 is supported");
        std::uint64_t total = gltfr_u32(data, 8);
        if (total > data.size())
            gltfr_fail("the GLB header declares more bytes than the file holds");
        std::size_t pos = 12;
        bool first = true;
        while (pos + 8 <= total) {
            const std::uint64_t len = gltfr_u32(data, pos);
            const std::uint32_t type = gltfr_u32(data, pos + 4);
            pos += 8;
            if (len > total - pos)
                gltfr_fail("a GLB chunk runs past the end of the file");
            if (first) {
                if (type != 0x4E4F534Au)
                    gltfr_fail("the first GLB chunk is not JSON");
                json_text.assign(data, pos, static_cast<std::size_t>(len));
                first = false;
            } else if (type == 0x004E4942u && !file.mHasBin) {
                file.mBin.assign(data, pos, static_cast<std::size_t>(len));
                file.mHasBin = true;
            }
            pos += static_cast<std::size_t>(len);
            if (len % 4 != 0)
                pos += static_cast<std::size_t>(4 - len % 4);
        }
        if (first)
            gltfr_fail("the GLB has no JSON chunk");
        while (!json_text.empty() && (json_text.back() == ' ' || json_text.back() == '\0'))
            json_text.pop_back();
    } else {
        if (gltfr_has_suffix_ci(rPath, ".glb"))
            gltfr_fail("'" + rPath + "' has the .glb suffix but not the GLB magic");
        json_text = data;
    }
    file.mRoot = detail::parse_json(json_text, "gltf");
    if (!file.mRoot.IsObject())
        gltfr_fail("the document is not a JSON object");
    const JsonValue* asset = gltfr_get(file.mRoot, "asset");
    const JsonValue* version = asset ? gltfr_get(*asset, "version") : nullptr;
    if (version == nullptr || !version->IsString() || version->mString.rfind("2.", 0) != 0)
        gltfr_fail("not a glTF 2.x document (asset.version)");
    return file;
}

std::string gltfr_percent_decode(const std::string& rText) {
    std::string out;
    for (std::size_t i = 0; i < rText.size(); ++i) {
        if (rText[i] == '%' && i + 2 < rText.size() && std::isxdigit(static_cast<unsigned char>(rText[i + 1])) &&
            std::isxdigit(static_cast<unsigned char>(rText[i + 2]))) {
            out += static_cast<char>(std::stoi(rText.substr(i + 1, 2), nullptr, 16));
            i += 2;
        } else {
            out += rText[i];
        }
    }
    return out;
}

// The bytes a buffer's `uri` names.
std::string gltfr_uri_bytes(const GltfrFile& rFile, const std::string& rUri) {
    if (rUri.rfind("data:", 0) == 0) {
        const std::size_t comma = rUri.find(',');
        const std::string header = comma == std::string::npos ? std::string() : rUri.substr(5, comma - 5);
        if (header.size() < 7 || header.compare(header.size() - 7, 7, ";base64") != 0)
            gltfr_fail("a data: URI that is not base64");
        const std::vector<unsigned char> bytes =
            detail::b64decode(rUri.data() + comma + 1, rUri.size() - comma - 1);
        return std::string(bytes.begin(), bytes.end());
    }
    const std::string path = gltfr_percent_decode(rUri);
    if (path.empty() || path[0] == '/' || path[0] == '\\' || path.find("://") != std::string::npos ||
        path.find('\\') != std::string::npos || path.find(':') != std::string::npos)
        gltfr_fail("the buffer URI '" + rUri + "' is not a relative file path");
    std::size_t from = 0;
    while (from <= path.size()) {
        std::size_t to = path.find('/', from);
        if (to == std::string::npos)
            to = path.size();
        if (path.compare(from, to - from, "..") == 0)
            gltfr_fail("the buffer URI '" + rUri + "' leaves the file's directory");
        from = to + 1;
    }
    return gltfr_read_file(rFile.mDirectory + path);
}

// ---- accessors ----------------------------------------------------------- //

struct GltfrBuffers {
    std::vector<std::string> mData;
};

GltfrBuffers gltfr_buffers(const GltfrFile& rFile) {
    GltfrBuffers out;
    const auto& list = gltfr_list(rFile.mRoot, "buffers");
    for (std::size_t i = 0; i < list.size(); ++i) {
        const JsonValue& b = gltfr_object(&list[i], "buffers[" + std::to_string(i) + "]");
        const std::uint64_t length = gltfr_uint_or(b, "byteLength", 0, "buffer");
        const JsonValue* uri = gltfr_get(b, "uri");
        std::string bytes;
        if (uri == nullptr) {
            if (i != 0 || !rFile.mHasBin)
                gltfr_fail("buffer " + std::to_string(i) + " has no uri and no GLB BIN chunk");
            bytes = rFile.mBin;
        } else {
            if (!uri->IsString())
                gltfr_fail("a buffer uri is not a string");
            bytes = gltfr_uri_bytes(rFile, uri->mString);
        }
        if (bytes.size() < length)
            gltfr_fail("buffer " + std::to_string(i) + " is shorter than its byteLength");
        out.mData.push_back(std::move(bytes));
    }
    return out;
}

int gltfr_component_size(std::uint64_t Type) {
    switch (Type) {
        case 5120:
        case 5121:
            return 1;
        case 5122:
        case 5123:
            return 2;
        case 5125:
        case 5126:
            return 4;
        default:
            gltfr_fail("unknown componentType " + std::to_string(Type));
    }
}

std::size_t gltfr_type_components(const std::string& rType) {
    if (rType == "SCALAR")
        return 1;
    if (rType == "VEC2")
        return 2;
    if (rType == "VEC3")
        return 3;
    if (rType == "VEC4" || rType == "MAT2")
        return 4;
    if (rType == "MAT3")
        return 9;
    if (rType == "MAT4")
        return 16;
    gltfr_fail("unknown accessor type '" + rType + "'");
}

double gltfr_component(const unsigned char* pBytes, std::uint64_t Type, bool Normalized) {
    switch (Type) {
        case 5120: {
            const double v = static_cast<std::int8_t>(pBytes[0]);
            return Normalized ? std::max(v / 127.0, -1.0) : v;
        }
        case 5121:
            return Normalized ? pBytes[0] / 255.0 : static_cast<double>(pBytes[0]);
        case 5122: {
            const double v = static_cast<std::int16_t>(pBytes[0] | (pBytes[1] << 8));
            return Normalized ? std::max(v / 32767.0, -1.0) : v;
        }
        case 5123: {
            const double v = static_cast<double>(pBytes[0] | (pBytes[1] << 8));
            return Normalized ? v / 65535.0 : v;
        }
        case 5125: {
            std::uint32_t u = 0;
            for (int k = 3; k >= 0; --k)
                u = (u << 8) | pBytes[k];
            return static_cast<double>(u);
        }
        default: {
            std::uint32_t u = 0;
            for (int k = 3; k >= 0; --k)
                u = (u << 8) | pBytes[k];
            float f;
            std::memcpy(&f, &u, sizeof(f));
            return static_cast<double>(f);
        }
    }
}

struct GltfrAccessor {
    std::size_t mCount = 0;
    std::size_t mNumComp = 0;
    std::string mType;
    std::uint64_t mComponentType = 0;
    bool mNormalized = false;
    std::string mName;
    std::vector<double> mValues;  // mCount * mNumComp
};

// The first byte of @p Count elements of @p ElemSize bytes at @p Stride in a view.
const unsigned char* gltfr_span(const GltfrFile& rFile, const GltfrBuffers& rBuffers,
                                std::uint64_t View, std::uint64_t ByteOffset, std::uint64_t ElemSize,
                                std::uint64_t Count, std::uint64_t& rStride, bool UseViewStride) {
    const auto& views = gltfr_list(rFile.mRoot, "bufferViews");
    if (View >= views.size())
        gltfr_fail("a bufferView index is out of range");
    const JsonValue& v = gltfr_object(&views[View], "bufferViews[" + std::to_string(View) + "]");
    if (const JsonValue* ext = gltfr_get(v, "extensions"))
        if (ext->Find("EXT_meshopt_compression") != nullptr)
            gltfr_fail("EXT_meshopt_compression is not supported");
    const std::uint64_t buffer = gltfr_uint_or(v, "buffer", std::uint64_t(-1), "bufferView");
    if (buffer >= rBuffers.mData.size())
        gltfr_fail("a bufferView names a buffer that does not exist");
    const std::uint64_t offset = gltfr_uint_or(v, "byteOffset", 0, "bufferView");
    const std::uint64_t length = gltfr_uint_or(v, "byteLength", 0, "bufferView");
    if (offset + length > rBuffers.mData[buffer].size())
        gltfr_fail("a bufferView runs past the end of its buffer");
    rStride = ElemSize;
    if (UseViewStride) {
        const std::uint64_t s = gltfr_uint_or(v, "byteStride", 0, "bufferView");
        if (s != 0) {
            if (s < ElemSize || s > (1u << 20))
                gltfr_fail("a bufferView byteStride is out of range");
            rStride = s;
        }
    }
    if (Count > 0) {
        const std::uint64_t last = ByteOffset + (Count - 1) * rStride + ElemSize;
        if (last > length)
            gltfr_fail("an accessor runs past the end of its bufferView");
    }
    return reinterpret_cast<const unsigned char*>(rBuffers.mData[buffer].data()) + offset +
           ByteOffset;
}

GltfrAccessor gltfr_read_accessor(const GltfrFile& rFile, const GltfrBuffers& rBuffers,
                                  std::uint64_t Index) {
    const auto& list = gltfr_list(rFile.mRoot, "accessors");
    if (Index >= list.size())
        gltfr_fail("an accessor index is out of range");
    const JsonValue& a = gltfr_object(&list[Index], "accessors[" + std::to_string(Index) + "]");
    GltfrAccessor out;
    const std::uint64_t count = gltfr_uint_or(a, "count", 0, "accessor");
    if (count > kGltfrMaxCount)
        gltfr_fail("an accessor count is too large");
    out.mCount = static_cast<std::size_t>(count);
    const JsonValue* type = gltfr_get(a, "type");
    if (type == nullptr || !type->IsString())
        gltfr_fail("an accessor has no type");
    out.mType = type->mString;
    out.mNumComp = gltfr_type_components(out.mType);
    out.mComponentType = gltfr_uint_or(a, "componentType", 0, "accessor");
    const int csize = gltfr_component_size(out.mComponentType);
    const JsonValue* norm = gltfr_get(a, "normalized");
    out.mNormalized = norm != nullptr && norm->IsBool() && norm->mBool;
    if (const JsonValue* name = gltfr_get(a, "name"))
        if (name->IsString())
            out.mName = name->mString;
    if (out.mType.rfind("MAT", 0) == 0)
        gltfr_fail("matrix accessors are not supported as vertex data");

    const std::uint64_t n = count * out.mNumComp;
    const JsonValue* view = gltfr_get(a, "bufferView");
    if (view == nullptr && n > kGltfrMaxZeros)
        gltfr_fail("an accessor without a bufferView is too large");
    out.mValues.assign(static_cast<std::size_t>(n), 0.0);
    if (view != nullptr) {
        const std::uint64_t elem = static_cast<std::uint64_t>(csize) * out.mNumComp;
        std::uint64_t stride = 0;
        const unsigned char* base =
            gltfr_span(rFile, rBuffers, gltfr_uint(*view, "accessor.bufferView"),
                       gltfr_uint_or(a, "byteOffset", 0, "accessor"), elem, count, stride, true);
        for (std::size_t i = 0; i < out.mCount; ++i)
            for (std::size_t c = 0; c < out.mNumComp; ++c)
                out.mValues[i * out.mNumComp + c] = gltfr_component(
                    base + i * stride + c * static_cast<std::size_t>(csize), out.mComponentType,
                    out.mNormalized);
    }
    if (const JsonValue* sparse = gltfr_get(a, "sparse")) {
        const std::uint64_t scount = gltfr_uint_or(*sparse, "count", 0, "sparse");
        if (scount > count)
            gltfr_fail("a sparse accessor has more entries than elements");
        const JsonValue& idx = gltfr_object(gltfr_get(*sparse, "indices"), "sparse.indices");
        const JsonValue& val = gltfr_object(gltfr_get(*sparse, "values"), "sparse.values");
        const std::uint64_t itype = gltfr_uint_or(idx, "componentType", 0, "sparse.indices");
        if (itype != 5121 && itype != 5123 && itype != 5125)
            gltfr_fail("sparse indices must be unsigned");
        const std::uint64_t isize = static_cast<std::uint64_t>(gltfr_component_size(itype));
        std::uint64_t s1 = 0, s2 = 0;
        const unsigned char* ib = gltfr_span(
            rFile, rBuffers, gltfr_uint_or(idx, "bufferView", std::uint64_t(-1), "sparse.indices"),
            gltfr_uint_or(idx, "byteOffset", 0, "sparse.indices"), isize, scount, s1, false);
        const std::uint64_t elem = static_cast<std::uint64_t>(csize) * out.mNumComp;
        const unsigned char* vb = gltfr_span(
            rFile, rBuffers, gltfr_uint_or(val, "bufferView", std::uint64_t(-1), "sparse.values"),
            gltfr_uint_or(val, "byteOffset", 0, "sparse.values"), elem, scount, s2, false);
        for (std::size_t i = 0; i < scount; ++i) {
            const double at = gltfr_component(ib + i * isize, itype, false);
            if (at >= static_cast<double>(count))
                gltfr_fail("a sparse index is out of range");
            for (std::size_t c = 0; c < out.mNumComp; ++c)
                out.mValues[static_cast<std::size_t>(at) * out.mNumComp + c] = gltfr_component(
                    vb + i * elem + c * static_cast<std::size_t>(csize), out.mComponentType,
                    out.mNormalized);
        }
    }
    return out;
}

// ---- transforms ---------------------------------------------------------- //

GltfrMat4 gltfr_identity() {
    GltfrMat4 m{};
    m[0] = m[5] = m[10] = m[15] = 1.0;
    return m;
}

GltfrMat4 gltfr_multiply(const GltfrMat4& rA, const GltfrMat4& rB) {
    GltfrMat4 m{};
    for (int col = 0; col < 4; ++col)
        for (int row = 0; row < 4; ++row) {
            double s = 0.0;
            for (int k = 0; k < 4; ++k)
                s += rA[static_cast<std::size_t>(k * 4 + row)] *
                     rB[static_cast<std::size_t>(col * 4 + k)];
            m[static_cast<std::size_t>(col * 4 + row)] = s;
        }
    return m;
}

GltfrMat4 gltfr_node_matrix(const JsonValue& rNode) {
    if (const JsonValue* m = gltfr_get(rNode, "matrix")) {
        const std::vector<double> v = gltfr_numbers(*m, 16, "node.matrix");
        GltfrMat4 out;
        std::copy(v.begin(), v.end(), out.begin());
        return out;
    }
    std::array<double, 3> t{0, 0, 0}, s{1, 1, 1};
    std::array<double, 4> q{0, 0, 0, 1};
    if (const JsonValue* v = gltfr_get(rNode, "translation")) {
        const auto n = gltfr_numbers(*v, 3, "node.translation");
        std::copy(n.begin(), n.end(), t.begin());
    }
    if (const JsonValue* v = gltfr_get(rNode, "scale")) {
        const auto n = gltfr_numbers(*v, 3, "node.scale");
        std::copy(n.begin(), n.end(), s.begin());
    }
    if (const JsonValue* v = gltfr_get(rNode, "rotation")) {
        const auto n = gltfr_numbers(*v, 4, "node.rotation");
        std::copy(n.begin(), n.end(), q.begin());
    }
    const double norm = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    if (!(norm > 0.0))
        gltfr_fail("a node rotation is the zero quaternion");
    const double x = q[0] / norm, y = q[1] / norm, z = q[2] / norm, w = q[3] / norm;
    GltfrMat4 m = gltfr_identity();
    m[0] = (1 - 2 * (y * y + z * z)) * s[0];
    m[1] = (2 * (x * y + z * w)) * s[0];
    m[2] = (2 * (x * z - y * w)) * s[0];
    m[4] = (2 * (x * y - z * w)) * s[1];
    m[5] = (1 - 2 * (x * x + z * z)) * s[1];
    m[6] = (2 * (y * z + x * w)) * s[1];
    m[8] = (2 * (x * z + y * w)) * s[2];
    m[9] = (2 * (y * z - x * w)) * s[2];
    m[10] = (1 - 2 * (x * x + y * y)) * s[2];
    m[12] = t[0];
    m[13] = t[1];
    m[14] = t[2];
    return m;
}

// The source-space offset of a root node `meshio++` as the writer emits it
// (an up-axis rotation, a uniform scale and a translation), or false.
bool gltfr_writer_root(const JsonValue& rNode, std::array<double, 3>& rOffset) {
    const JsonValue* name = gltfr_get(rNode, "name");
    if (name == nullptr || !name->IsString() || name->mString != "meshio++" ||
        gltfr_get(rNode, "mesh") != nullptr || gltfr_get(rNode, "matrix") != nullptr)
        return false;
    std::array<double, 3> t{0, 0, 0};
    double s = 1.0;
    int up = 2;  // 0 = Z up, 1 = X up, 2 = Y up (no rotation)
    if (const JsonValue* v = gltfr_get(rNode, "rotation")) {
        const auto q = gltfr_numbers(*v, 4, "node.rotation");
        auto close = [&](double a, double b, double c, double d) {
            const double tol = 1e-9;
            const bool same = std::abs(q[0] - a) < tol && std::abs(q[1] - b) < tol &&
                              std::abs(q[2] - c) < tol && std::abs(q[3] - d) < tol;
            const bool flip = std::abs(q[0] + a) < tol && std::abs(q[1] + b) < tol &&
                              std::abs(q[2] + c) < tol && std::abs(q[3] + d) < tol;
            return same || flip;
        };
        if (close(-kGltfrHalfSqrt2, 0, 0, kGltfrHalfSqrt2))
            up = 0;
        else if (close(0, 0, kGltfrHalfSqrt2, kGltfrHalfSqrt2))
            up = 1;
        else if (!close(0, 0, 0, 1))
            return false;
    }
    if (const JsonValue* v = gltfr_get(rNode, "scale")) {
        const auto n = gltfr_numbers(*v, 3, "node.scale");
        if (!(n[0] > 0.0) || n[1] != n[0] || n[2] != n[0])
            return false;
        s = n[0];
    }
    if (const JsonValue* v = gltfr_get(rNode, "translation")) {
        const auto n = gltfr_numbers(*v, 3, "node.translation");
        std::copy(n.begin(), n.end(), t.begin());
    }
    // center = R^-1 t / s (the writer's translation is R * (s * center)).
    if (up == 0)
        rOffset = {t[0] / s, -t[2] / s, t[1] / s};
    else if (up == 1)
        rOffset = {t[1] / s, -t[0] / s, t[2] / s};
    else
        rOffset = {t[0] / s, t[1] / s, t[2] / s};
    for (double& o : rOffset)
        o += 0.0;
    return true;
}

// ---- the mesh ------------------------------------------------------------ //

enum GltfrKind { kGltfrTri = 0, kGltfrLine = 1, kGltfrVertex = 2 };

struct GltfrField {
    std::size_t mNumComp = 0;
    std::vector<double> mData;
    bool mDropped = false;
};

struct GltfrState {
    std::vector<double> mPoints;  // 3 per vertex
    std::vector<std::int64_t> mCells[3];
    std::map<std::string, std::array<std::vector<std::int64_t>, 3>> mRegions;
    std::map<std::string, GltfrField> mFields;
    std::set<std::string> mWarned;
    bool mWantData = true;
    const ReadOptions* mpRead = nullptr;
};

std::string gltfr_field_name(const std::string& rAttr, const GltfrAccessor& rAcc, int Mode) {
    if (rAttr == "COLOR_0")
        return "gltf:color";
    if (rAttr == "NORMAL")
        return Mode == 0 ? "normals" : std::string();
    if (rAttr.size() > 1 && rAttr[0] == '_')
        return rAcc.mName.empty() ? rAttr.substr(1) : rAcc.mName;
    return std::string();
}

void gltfr_primitive(const GltfrFile& rFile, const GltfrBuffers& rBuffers, const JsonValue& rPrim,
                     const GltfrMat4& rWorld, const std::array<double, 3>& rOffset,
                     const std::string& rRegion, GltfrState& rState) {
    if (const JsonValue* ext = gltfr_get(rPrim, "extensions"))
        if (ext->Find("KHR_draco_mesh_compression") != nullptr)
            gltfr_fail("KHR_draco_mesh_compression is not supported");
    const std::uint64_t mode = gltfr_uint_or(rPrim, "mode", 4, "primitive");
    if (mode > 6)
        gltfr_fail("unknown primitive mode " + std::to_string(mode));
    const JsonValue& attrs = gltfr_object(gltfr_get(rPrim, "attributes"), "primitive.attributes");
    const JsonValue* position = attrs.Find("POSITION");
    if (position == nullptr) {
        if (rState.mWarned.insert("noposition").second)
            log::warn("gltf: a primitive without POSITION is skipped.");
        return;
    }
    const GltfrAccessor pos =
        gltfr_read_accessor(rFile, rBuffers, gltfr_uint(*position, "POSITION"));
    if (pos.mNumComp != 3)
        gltfr_fail("POSITION is not VEC3");
    const std::size_t nv = pos.mCount;
    const std::size_t base = rState.mPoints.size() / 3;
    if (base + nv > kGltfrMaxCount)
        gltfr_fail("the file holds too many vertices");

    // indices
    std::vector<std::uint64_t> index;
    if (const JsonValue* ii = gltfr_get(rPrim, "indices")) {
        const GltfrAccessor acc = gltfr_read_accessor(rFile, rBuffers, gltfr_uint(*ii, "indices"));
        if (acc.mNumComp != 1 || (acc.mComponentType != 5121 && acc.mComponentType != 5123 &&
                                  acc.mComponentType != 5125))
            gltfr_fail("indices must be unsigned SCALAR");
        index.reserve(acc.mCount);
        for (double v : acc.mValues) {
            if (v >= static_cast<double>(nv))
                gltfr_fail("an index is out of range");
            index.push_back(static_cast<std::uint64_t>(v));
        }
    } else {
        index.resize(nv);
        for (std::size_t i = 0; i < nv; ++i)
            index[i] = i;
    }

    // vertices, in world space
    for (std::size_t i = 0; i < nv; ++i) {
        const double x = pos.mValues[i * 3], y = pos.mValues[i * 3 + 1], z = pos.mValues[i * 3 + 2];
        for (int r = 0; r < 3; ++r) {
            const double w = rWorld[static_cast<std::size_t>(r)] * x +
                             rWorld[static_cast<std::size_t>(4 + r)] * y +
                             rWorld[static_cast<std::size_t>(8 + r)] * z +
                             rWorld[static_cast<std::size_t>(12 + r)] + rOffset[static_cast<std::size_t>(r)];
            if (!std::isfinite(w))
                gltfr_fail("a vertex position is not finite");
            rState.mPoints.push_back(w + 0.0);
        }
    }

    // attributes
    std::map<std::string, GltfrField> local;
    if (rState.mWantData) {
        for (const auto& member : attrs.mMembers) {
            if (member.mKey == "POSITION")
                continue;
            if (!member.mValue.IsNumber())
                continue;
            if (member.mKey.empty() || (member.mKey[0] != '_' && member.mKey != "COLOR_0" &&
                                        member.mKey != "NORMAL"))
                continue;
            const GltfrAccessor acc =
                gltfr_read_accessor(rFile, rBuffers, gltfr_uint(member.mValue, member.mKey));
            if (acc.mCount != nv)
                gltfr_fail("attribute " + member.mKey + " does not have one value per vertex");
            const std::string name = gltfr_field_name(member.mKey, acc, static_cast<int>(mode));
            if (name.empty() || acc.mNumComp > 4 || !rState.mpRead->WantsArray(name))
                continue;
            if (!local.emplace(name, GltfrField{acc.mNumComp, acc.mValues, false}).second &&
                rState.mWarned.insert("dup:" + name).second)
                log::warn("gltf: two attributes map to '{}'; the first is kept.", name);
        }
    }
    for (auto& [name, field] : local) {
        auto it = rState.mFields.find(name);
        if (it == rState.mFields.end()) {
            GltfrField f;
            f.mNumComp = field.mNumComp;
            f.mData.assign(base * field.mNumComp, std::nan(""));
            it = rState.mFields.emplace(name, std::move(f)).first;
        } else if (it->second.mNumComp != field.mNumComp && !it->second.mDropped) {
            it->second.mDropped = true;
            log::warn(
                "gltf: attribute '{}' changes its component count between primitives; dropped.",
                name);
        }
    }
    for (auto& [name, f] : rState.mFields) {
        f.mData.resize((base + nv) * f.mNumComp, std::nan(""));
        auto it = local.find(name);
        if (it != local.end() && it->second.mNumComp == f.mNumComp)
            std::copy(it->second.mData.begin(), it->second.mData.end(),
                      f.mData.begin() + static_cast<std::ptrdiff_t>(base * f.mNumComp));
    }

    // cells
    auto add = [&](GltfrKind kind, std::initializer_list<std::uint64_t> verts) {
        for (std::uint64_t v : verts)
            rState.mCells[kind].push_back(static_cast<std::int64_t>(base + v));
    };
    const std::size_t before[3] = {rState.mCells[0].size() / 3, rState.mCells[1].size() / 2,
                                   rState.mCells[2].size()};
    const std::size_t n = index.size();
    switch (mode) {
        case 0:
            for (std::size_t i = 0; i < n; ++i)
                add(kGltfrVertex, {index[i]});
            break;
        case 1:
            if (n % 2 != 0)
                gltfr_fail("a LINES primitive has an odd number of indices");
            for (std::size_t i = 0; i + 1 < n; i += 2)
                add(kGltfrLine, {index[i], index[i + 1]});
            break;
        case 2:
        case 3:
            for (std::size_t i = 0; i + 1 < n; ++i)
                add(kGltfrLine, {index[i], index[i + 1]});
            if (mode == 2 && n > 2)
                add(kGltfrLine, {index[n - 1], index[0]});
            break;
        case 4:
            if (n % 3 != 0)
                gltfr_fail("a TRIANGLES primitive has a count that is not a multiple of 3");
            for (std::size_t i = 0; i + 2 < n; i += 3)
                add(kGltfrTri, {index[i], index[i + 1], index[i + 2]});
            break;
        case 5:
            for (std::size_t i = 0; i + 2 < n; ++i) {
                const std::uint64_t a = index[i], b = index[i + 1], c = index[i + 2];
                if (a == b || b == c || a == c)
                    continue;  // degenerate joins
                if (i % 2 == 0)
                    add(kGltfrTri, {a, b, c});
                else
                    add(kGltfrTri, {b, a, c});
            }
            break;
        default:
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::uint64_t a = index[0], b = index[i], c = index[i + 1];
                if (a == b || b == c || a == c)
                    continue;
                add(kGltfrTri, {a, b, c});
            }
            break;
    }
    if (!rRegion.empty() && rRegion != "unassigned" && rRegion != "mesh") {
        auto& entry = rState.mRegions[rRegion];
        const std::size_t after[3] = {rState.mCells[0].size() / 3, rState.mCells[1].size() / 2,
                                      rState.mCells[2].size()};
        for (int k = 0; k < 3; ++k)
            for (std::size_t c = before[k]; c < after[k]; ++c)
                entry[static_cast<std::size_t>(k)].push_back(static_cast<std::int64_t>(c));
    }
}

void gltfr_walk(const GltfrFile& rFile, const GltfrBuffers& rBuffers, std::size_t NodeIndex,
                const GltfrMat4& rParent, const std::array<double, 3>& rOffset,
                std::vector<char>& rVisited, GltfrState& rState) {
    const auto& nodes = gltfr_list(rFile.mRoot, "nodes");
    if (NodeIndex >= nodes.size())
        gltfr_fail("a node index is out of range");
    if (rVisited[NodeIndex])
        gltfr_fail("the node graph is not a forest (a node is reached twice)");
    rVisited[NodeIndex] = 1;
    const JsonValue& node = gltfr_object(&nodes[NodeIndex], "nodes[" + std::to_string(NodeIndex) + "]");
    const GltfrMat4 world = gltfr_multiply(rParent, gltfr_node_matrix(node));
    if (const JsonValue* mesh = gltfr_get(node, "mesh")) {
        const auto& meshes = gltfr_list(rFile.mRoot, "meshes");
        const std::uint64_t mi = gltfr_uint(*mesh, "node.mesh");
        if (mi >= meshes.size())
            gltfr_fail("a mesh index is out of range");
        const JsonValue& m = gltfr_object(&meshes[mi], "meshes[" + std::to_string(mi) + "]");
        const JsonValue* name = gltfr_get(node, "name");
        const std::string region = name && name->IsString() ? name->mString : std::string();
        for (const JsonValue& prim : gltfr_list(m, "primitives")) {
            if (gltfr_get(prim, "targets") != nullptr && rState.mWarned.insert("morph").second)
                log::warn("gltf: morph targets are skipped.");
            gltfr_primitive(rFile, rBuffers, gltfr_object(&prim, "primitive"), world, rOffset,
                            region, rState);
        }
    }
    if (const JsonValue* children = gltfr_get(node, "children")) {
        if (!children->IsArray())
            gltfr_fail("node.children is not an array");
        for (const JsonValue& c : children->mItems)
            gltfr_walk(rFile, rBuffers, static_cast<std::size_t>(gltfr_uint(c, "node.children")),
                       world, rOffset, rVisited, rState);
    }
}

}  // namespace

Mesh read_gltf(const std::string& rPath, const GltfReadOptions& rOptions,
               const ReadOptions& rReadOptions) {
    const GltfrFile file = gltfr_load(rPath);
    for (const JsonValue& ext : gltfr_list(file.mRoot, "extensionsRequired")) {
        if (!ext.IsString())
            continue;
        if (ext.mString == "KHR_draco_mesh_compression" || ext.mString == "EXT_meshopt_compression")
            gltfr_fail(ext.mString + " is not supported");
    }
    const GltfrBuffers buffers = gltfr_buffers(file);

    GltfrState state;
    state.mWantData = rReadOptions.WantsAnyData();
    state.mpRead = &rReadOptions;

    std::vector<std::string> skipped;
    for (const char* key : {"animations", "skins"})
        if (!gltfr_list(file.mRoot, key).empty())
            skipped.emplace_back(key);
    if (!gltfr_list(file.mRoot, "textures").empty() || !gltfr_list(file.mRoot, "images").empty())
        skipped.emplace_back("textures");
    if (!skipped.empty()) {
        std::string joined;
        for (std::size_t i = 0; i < skipped.size(); ++i)
            joined += (i ? ", " : "") + skipped[i];
        log::warn("gltf: {} are skipped.", joined);
    }

    // The scene's roots.
    const auto& nodes = gltfr_list(file.mRoot, "nodes");
    std::vector<std::size_t> roots;
    const auto& scenes = gltfr_list(file.mRoot, "scenes");
    if (!scenes.empty()) {
        const std::uint64_t si = gltfr_uint_or(file.mRoot, "scene", 0, "document");
        if (si >= scenes.size())
            gltfr_fail("the default scene index is out of range");
        for (const JsonValue& n : gltfr_list(gltfr_object(&scenes[si], "scene"), "nodes"))
            roots.push_back(static_cast<std::size_t>(gltfr_uint(n, "scene.nodes")));
    } else {
        std::vector<char> is_child(nodes.size(), 0);
        for (const JsonValue& n : nodes)
            for (const JsonValue& c : gltfr_list(n, "children")) {
                const std::uint64_t ci = gltfr_uint(c, "node.children");
                if (ci < nodes.size())
                    is_child[static_cast<std::size_t>(ci)] = 1;
            }
        for (std::size_t i = 0; i < nodes.size(); ++i)
            if (!is_child[i])
                roots.push_back(i);
    }

    std::vector<char> visited(nodes.size(), 0);
    for (std::size_t root : roots) {
        if (root >= nodes.size())
            gltfr_fail("a scene names a node that does not exist");
        std::array<double, 3> offset{0, 0, 0};
        const JsonValue& rootNode = gltfr_object(&nodes[root], "node");
        if (roots.size() == 1 && gltfr_writer_root(rootNode, offset)) {
            visited[root] = 1;
            GltfrMat4 shift = gltfr_identity();
            for (const JsonValue& c : gltfr_list(rootNode, "children"))
                gltfr_walk(file, buffers, static_cast<std::size_t>(gltfr_uint(c, "node.children")),
                           shift, offset, visited, state);
        } else {
            gltfr_walk(file, buffers, root, gltfr_identity(), offset, visited, state);
        }
    }

    // Welding.
    std::size_t n = state.mPoints.size() / 3;
    if (rOptions.mWeld && n > 0) {
        struct Key {
            std::array<double, 3> v;
            bool operator==(const Key& o) const { return v == o.v; }
        };
        struct KeyHash {
            std::size_t operator()(const Key& k) const {
                std::size_t h = 1469598103934665603ull;
                for (double d : k.v) {
                    std::uint64_t b;
                    std::memcpy(&b, &d, sizeof(b));
                    h = (h ^ b) * 1099511628211ull;
                }
                return h;
            }
        };
        std::unordered_map<Key, std::int64_t, KeyHash> seen;
        std::vector<std::int64_t> remap(n);
        std::vector<std::size_t> first;
        for (std::size_t i = 0; i < n; ++i) {
            const Key k{{state.mPoints[i * 3], state.mPoints[i * 3 + 1], state.mPoints[i * 3 + 2]}};
            const auto it = seen.emplace(k, static_cast<std::int64_t>(first.size()));
            if (it.second)
                first.push_back(i);
            remap[i] = it.first->second;
        }
        std::vector<double> points;
        points.reserve(first.size() * 3);
        for (std::size_t i : first)
            points.insert(points.end(), state.mPoints.begin() + static_cast<std::ptrdiff_t>(i * 3),
                          state.mPoints.begin() + static_cast<std::ptrdiff_t>(i * 3 + 3));
        for (auto& [name, f] : state.mFields) {
            std::vector<double> data;
            data.reserve(first.size() * f.mNumComp);
            for (std::size_t i : first)
                data.insert(data.end(), f.mData.begin() + static_cast<std::ptrdiff_t>(i * f.mNumComp),
                            f.mData.begin() + static_cast<std::ptrdiff_t>((i + 1) * f.mNumComp));
            f.mData = std::move(data);
        }
        for (auto& cells : state.mCells)
            for (std::int64_t& v : cells)
                v = remap[static_cast<std::size_t>(v)];
        state.mPoints = std::move(points);
        n = first.size();
    }

    Mesh mesh;
    NDArray pts(DType::Float64, {n, static_cast<std::size_t>(3)});
    std::copy(state.mPoints.begin(), state.mPoints.end(), pts.As<double>());
    mesh.AssignPoints(std::move(pts));

    static const char* const kTypes[3] = {"triangle", "line", "vertex"};
    static const std::size_t kWidth[3] = {3, 2, 1};
    std::int64_t offsets[3] = {0, 0, 0};
    std::int64_t total = 0;
    for (int k = 0; k < 3; ++k) {
        const std::size_t cells = state.mCells[k].size() / kWidth[k];
        offsets[k] = total;
        if (cells == 0)
            continue;
        total += static_cast<std::int64_t>(cells);
        NDArray conn(DType::Int64, {cells, kWidth[k]});
        std::copy(state.mCells[k].begin(), state.mCells[k].end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(kTypes[k], std::move(conn));
    }

    for (auto& [name, f] : state.mFields) {
        if (f.mDropped)
            continue;
        NDArray arr = f.mNumComp == 1 ? NDArray(DType::Float64, {n})
                                      : NDArray(DType::Float64, {n, f.mNumComp});
        std::copy(f.mData.begin(), f.mData.end(), arr.As<double>());
        mesh.AddPointData(name, std::move(arr));
    }

    for (auto& [name, entry] : state.mRegions) {
        std::vector<std::int64_t> ids;
        for (int k = 0; k < 3; ++k)
            for (std::int64_t c : entry[static_cast<std::size_t>(k)])
                ids.push_back(offsets[k] + c);
        NDArray arr(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), arr.As<std::int64_t>());
        mesh.AddRegion(Region(name, RegionKind::Cell, std::move(arr)));
    }
    return mesh;
}

}  // namespace meshioplusplus
