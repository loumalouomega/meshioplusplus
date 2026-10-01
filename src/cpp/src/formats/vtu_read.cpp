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
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "../detail/region_field_data.hpp"
#include "vtk_preflight.hpp"
#include "../detail/vtu_decode.hpp"
#include "../detail/open_source.hpp"
#include "../detail/vtk_xml_read.hpp"

namespace meshioplusplus {

namespace detail {
namespace {

using detail::vtu_to_int64;

/**
 * @brief Sequential reader over one array's encoded bytes: either raw bytes (a
 * raw `<AppendedData>` payload) or base64 text.
 *
 * Base64 is decoded group by group, so a header and a body that were encoded
 * separately (each padded, as VTK's appended writer does) and one stream that
 * encodes both read the same.
 */
struct VtuByteSource {
    const unsigned char* mRaw = nullptr;
    std::size_t mRawLen = 0;
    const char* mText = nullptr;
    std::size_t mTextLen = 0;
    std::size_t mPos = 0;
    std::vector<unsigned char> mPending;
    std::size_t mPendingPos = 0;

    static int Base64Value(char c) {
        if (c >= 'A' && c <= 'Z')
            return c - 'A';
        if (c >= 'a' && c <= 'z')
            return c - 'a' + 26;
        if (c >= '0' && c <= '9')
            return c - '0' + 52;
        if (c == '+')
            return 62;
        if (c == '/')
            return 63;
        return -1;
    }

    // Upper bound on the bytes left, to refuse absurd sizes before allocating.
    std::size_t Remaining() const {
        if (mRaw)
            return mRawLen - std::min(mPos, mRawLen);
        return (mPending.size() - mPendingPos) + (mTextLen - std::min(mPos, mTextLen)) / 4 * 3 + 3;
    }

    void Take(unsigned char* pOut, std::size_t n) {
        if (mRaw) {
            if (mPos > mRawLen || n > mRawLen - mPos)
                throw ReadError("VTU: appended array runs past the end of the data");
            std::memcpy(pOut, mRaw + mPos, n);
            mPos += n;
            return;
        }
        mPending.reserve(mPendingPos + n + 3);
        while (mPending.size() - mPendingPos < n) {
            char q[4];
            int k = 0;
            while (k < 4 && mPos < mTextLen) {
                const char c = mText[mPos++];
                if (!std::isspace(static_cast<unsigned char>(c)))
                    q[k++] = c;
            }
            if (k < 4)
                throw ReadError("VTU: base64 data ends early");
            std::uint32_t t = 0;
            int pad = 0;
            for (int i = 0; i < 4; ++i) {
                int v = 0;
                if (q[i] == '=')
                    ++pad;
                else if ((v = Base64Value(q[i])) < 0)
                    throw ReadError("VTU: invalid base64 character");
                t = (t << 6) | static_cast<std::uint32_t>(v);
            }
            mPending.push_back(static_cast<unsigned char>((t >> 16) & 0xff));
            if (pad < 2)
                mPending.push_back(static_cast<unsigned char>((t >> 8) & 0xff));
            if (pad < 1)
                mPending.push_back(static_cast<unsigned char>(t & 0xff));
        }
        std::memcpy(pOut, mPending.data() + mPendingPos, n);
        mPendingPos += n;
        if (mPendingPos == mPending.size()) {
            mPending.clear();
            mPendingPos = 0;
        }
    }
};

std::uint64_t vtu_read_uint(const unsigned char* pP, std::size_t Size, bool BigEndian) {
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < Size; ++i)
        v = BigEndian ? (v << 8) | pP[i] : v | (static_cast<std::uint64_t>(pP[i]) << (8 * i));
    return v;
}

/** @brief One array's bytes: the byte-count header (or, compressed, the block
 * table) and the body, decompressed; the header is in the file's byte order. */
std::vector<unsigned char> vtu_decode_sequential(VtuByteSource& rSrc, detail::VtkCodec Codec,
                                                 std::size_t HeaderSize, bool BigEndian) {
    unsigned char h[8];
    auto next = [&]() {
        rSrc.Take(h, HeaderSize);
        return vtu_read_uint(h, HeaderSize, BigEndian);
    };
    std::vector<unsigned char> out;
    if (Codec == detail::VtkCodec::None) {
        const std::uint64_t n = next();
        if (n > rSrc.Remaining())
            throw ReadError("VTU: array size exceeds the data");
        out.resize(static_cast<std::size_t>(n));
        if (n)
            rSrc.Take(out.data(), out.size());
        return out;
    }
    const std::uint64_t num_blocks = next();
    const std::uint64_t max_block = next();
    const std::uint64_t last_block = next();
    if (num_blocks > rSrc.Remaining() / HeaderSize)
        throw ReadError("VTU: block count exceeds the data");
    std::vector<std::uint64_t> sizes(static_cast<std::size_t>(num_blocks));
    for (auto& rSize : sizes)
        rSize = next();
    // The blocks' decompressed total, reserved when it is plausible: no codec
    // expands a block by more than about 2^16, as vtk_codec_decompress_block
    // enforces per block.
    std::uint64_t comp_total = 0;
    for (const std::uint64_t z : sizes) {
        // The compressed blocks are read from the data that follows, so their
        // total cannot exceed it.
        if (z > rSrc.Remaining() || comp_total > rSrc.Remaining() - z)
            throw ReadError("VTU: compressed blocks exceed the data");
        comp_total += z;
    }
    const std::uint64_t ceiling = (std::uint64_t{1} << 16) * (comp_total + 1) + (1u << 20);
    if (num_blocks > 0 && max_block <= ceiling && last_block <= ceiling &&
        num_blocks - 1 <= ceiling / std::max<std::uint64_t>(max_block, 1)) {
        const std::uint64_t want = (num_blocks - 1) * max_block + last_block;
        if (want <= ceiling)
            out.reserve(static_cast<std::size_t>(std::min<std::uint64_t>(want, 1u << 26)));
    }
    std::vector<unsigned char> comp;
    for (std::size_t k = 0; k < sizes.size(); ++k) {
        if (sizes[k] > rSrc.Remaining())
            throw ReadError("VTU: compressed block exceeds the data");
        comp.resize(static_cast<std::size_t>(sizes[k]));
        if (!comp.empty())
            rSrc.Take(comp.data(), comp.size());
        const std::size_t expected =
            static_cast<std::size_t>(k + 1 == sizes.size() ? last_block : max_block);
        std::vector<unsigned char> dec =
            detail::vtk_codec_decompress_block(Codec, comp.data(), comp.size(), expected);
        if (dec.size() != expected)
            throw ReadError("VTK XML: decompressed block size differs from its header");
        out.insert(out.end(), dec.begin(), dec.end());
    }
    return out;
}

NDArray vtu_array_from_bytes(const std::vector<unsigned char>& rBytes, DType Dt, bool BigEndian) {
    const std::size_t isz = dtype_size(Dt);
    if (!isz || rBytes.size() % isz != 0)
        throw ReadError("VTK XML: array byte count is not a multiple of its element size");
    const std::size_t n = isz ? rBytes.size() / isz : 0;
    NDArray a(Dt, {n});
    if (n)
        std::memcpy(a.Data(), rBytes.data(), n * isz);
    if (BigEndian && isz > 1) {
        char* p = reinterpret_cast<char*>(a.Data());
        for (std::size_t i = 0; i < n; ++i)
            detail::bswap_inplace(p + i * isz, static_cast<int>(isz));
    }
    return a;
}

}  // namespace

NDArray vtu_read_data_array(const pugi::xml_node& rDa, const VtuContext& rCtx,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);
    if (rDa.attribute("NumberOfComponents") && rNumComponents <= 0)
        throw ReadError("VTK XML: NumberOfComponents must be positive");

    auto decode = [&]() -> NDArray {
        if (fmt == "ascii")
            return detail::vtu_parse_ascii(rDa.text().get(), dt);
        if (fmt == "binary") {
            if (!rCtx.mBigEndian)
                return detail::vtu_decode_bin_view(detail::vtu_strip_view(rDa.text().get()), dt,
                                                   rCtx.mCodec, rCtx.mHeaderSize);
            const char* text = rDa.text().get();
            VtuByteSource src;
            src.mText = text;
            src.mTextLen = std::strlen(text);
            return vtu_array_from_bytes(
                vtu_decode_sequential(src, rCtx.mCodec, rCtx.mHeaderSize, true), dt, true);
        }
        if (fmt == "appended") {
            // Accept padding, but not missing, negative, overflowing or junk offsets.
            const char* text = rDa.attribute("offset").as_string("");
            while (*text && std::isspace(static_cast<unsigned char>(*text)))
                ++text;
            if (*text < '0' || *text > '9')
                throw ReadError("VTK XML: invalid appended offset");
            char* end = nullptr;
            errno = 0;
            const std::uint64_t offset = std::strtoull(text, &end, 10);
            if (errno == ERANGE)
                throw ReadError("VTK XML: appended offset overflows UInt64");
            while (*end && std::isspace(static_cast<unsigned char>(*end)))
                ++end;
            if (*end)
                throw ReadError("VTK XML: invalid appended offset");
            VtuByteSource src;
            if (rCtx.mRaw) {
                src.mRaw = rCtx.mRaw;
                src.mRawLen = rCtx.mRawLen;
            } else if (rCtx.mBase64) {
                src.mText = rCtx.mBase64;
                src.mTextLen = rCtx.mBase64Len;
            } else {
                throw ReadError("VTU: appended DataArray but no <AppendedData>");
            }
            if (offset > (src.mRaw ? src.mRawLen : src.mTextLen))
                throw ReadError("VTU: appended offset past the end of the data");
            src.mPos = static_cast<std::size_t>(offset);
            return vtu_array_from_bytes(
                vtu_decode_sequential(src, rCtx.mCodec, rCtx.mHeaderSize, rCtx.mBigEndian), dt,
                rCtx.mBigEndian);
        }
        throw ReadError("VTU '" + fmt + "' data is not supported by the C++ reader");
    };
    NDArray out = decode();
    if (rNumComponents > 0 && out.Size() % static_cast<std::size_t>(rNumComponents) != 0)
        throw ReadError("VTK XML: array size does not fit NumberOfComponents");
    return out;
}

/**
 * @brief The XML document and the file's bytes when they are needed.
 *
 * A raw `<AppendedData>` payload is not XML text (it may hold any byte, `<`
 * included): the file is read as bytes, the payload cut out, and only the text
 * around it parsed. Every other file parses as before.
 */
void vtu_load(const std::string& rPath, unsigned int ParseOptions, VtuSource& rSource,
              const char* pType, const char* pFormat) {
    const std::string format(pFormat);
    const std::string decline = "lzma-compressed " + format + " not supported by the C++ reader";
    detail::vtk_preflight(rPath, pType, decline.c_str());
    // The file is read once: a raw <AppendedData> payload is not XML, so only
    // the text around it is parsed, and the payload is a view into the same
    // bytes (it read the file twice and parsed it twice before v16.21.0).
    rSource.mFile.emplace(detail::open_source(rPath, "Could not open file: " + rPath));
    rSource.mBytes = rSource.mFile->View();
    const std::string_view b = rSource.mBytes;
    const std::size_t tag = b.find("<AppendedData");
    const std::size_t tag_end = tag == std::string_view::npos ? tag : b.find('>', tag);
    const bool raw =
        tag_end != std::string_view::npos &&
        detail::vtk_preflight_attribute(b.substr(tag, tag_end - tag), "encoding") == "raw";
    pugi::xml_parse_result res;
    if (!raw) {
        res = rSource.mDoc.load_buffer(b.data(), b.size(), ParseOptions);
        if (!res)
            throw ReadError(format + " XML parse failed: " + res.description());
        return;
    }
    std::size_t underscore = tag_end + 1;
    while (underscore < b.size() && std::isspace(static_cast<unsigned char>(b[underscore])))
        ++underscore;
    const std::size_t stop = b.rfind("</AppendedData>");
    if (underscore >= b.size() || b[underscore] != '_' || stop == std::string::npos ||
        stop <= underscore)
        throw ReadError(format + ": AppendedData must start with '_' and be closed");
    std::string xml(b.substr(0, tag_end + 1));
    xml += b.substr(stop);
    rSource.mDoc.reset();
    res = rSource.mDoc.load_buffer(xml.data(), xml.size(), ParseOptions);
    if (!res)
        throw ReadError(format + " XML parse failed: " + res.description());
    rSource.mIsRaw = true;
    rSource.mRawStart = underscore + 1;
    rSource.mRawStop = stop;
}

VtuContext vtk_xml_read_context(const VtuSource& rSource, const char* pFormat) {
    VtuContext ctx;
    const std::string format(pFormat);
    const auto root = rSource.mDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    const std::string compressor = root.attribute("compressor").as_string("");
    if (!compressor.empty()) {
        bool found = false;
        for (const auto codec : {VtkCodec::Zlib, VtkCodec::LZ4, VtkCodec::ZSTD})
            if (compressor == vtk_codec_compressor(codec)) {
                ctx.mCodec = codec;
                found = true;
            }
        if (!found) {
            if (compressor == vtk_codec_compressor(VtkCodec::LZMA))
                throw ReadError("lzma-compressed " + format + " not supported by the C++ reader");
            throw ReadError("Unknown " + format + " compressor '" + compressor + "'");
        }
    }
    vtk_codec_require_read(ctx.mCodec);
    const std::string header = root.attribute("header_type").as_string("UInt32");
    if (header != "UInt32" && header != "UInt64")
        throw ReadError("Unknown " + format + " header type '" + header + "'");
    ctx.mHeaderSize = header == "UInt64" ? 8 : 4;
    const std::string order = root.attribute("byte_order").as_string("LittleEndian");
    if (order != "LittleEndian" && order != "BigEndian")
        throw ReadError("Unknown " + format + " byte order '" + order + "'");
    ctx.mBigEndian = order == "BigEndian";
    if (rSource.mIsRaw) {
        ctx.mRaw =
            reinterpret_cast<const unsigned char*>(rSource.mBytes.data()) + rSource.mRawStart;
        ctx.mRawLen = rSource.mRawStop - rSource.mRawStart;
    } else if (const auto app = root.child("AppendedData")) {
        const std::string encoding = app.attribute("encoding").as_string("base64");
        if (encoding != "base64")
            throw ReadError("Unknown " + format + " AppendedData encoding '" + encoding + "'");
        const char* text = app.text().get();
        while (*text && std::isspace(static_cast<unsigned char>(*text)))
            ++text;
        if (*text != '_')
            throw ReadError(format + ": AppendedData does not start with '_'");
        ctx.mBase64 = text + 1;
        ctx.mBase64Len = std::strlen(text + 1);
    }
    return ctx;
}

}  // namespace detail

namespace {
using detail::vtu_load;
using detail::vtu_read_data_array;
using detail::vtu_to_int64;
using detail::VtuContext;
using detail::VtuSource;

/**
 * @brief The `<Piece>` nodes plus the framing attributes every path needs.
 *
 * Shared by the mesh and metadata readers so the two cannot disagree about
 * which files they accept -- a metadata summary must never succeed on a file
 * `read_vtu` would reject.
 */
struct vtu_header {
    pugi::xml_node mGrid;
    std::vector<pugi::xml_node> mPieces;
    VtuContext mCtx;
    std::size_t mNumPoints = 0;
};

vtu_header vtu_parse_header(const VtuSource& rSource) {
    pugi::xml_node root = rSource.mDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "UnstructuredGrid")
        throw ReadError("Expected type UnstructuredGrid");

    vtu_header h;
    h.mCtx = detail::vtk_xml_read_context(rSource, "VTU");

    pugi::xml_node grid = root.child("UnstructuredGrid");
    if (!grid)
        throw ReadError("No UnstructuredGrid found");

    h.mGrid = grid;
    for (pugi::xml_node piece : grid.children("Piece")) {
        const std::size_t n =
            static_cast<std::size_t>(piece.attribute("NumberOfPoints").as_ullong());
        if (n > 0 && !piece.child("Points").child("DataArray"))
            throw ReadError("VTU: a Piece declares points but has no <Points>");
        h.mPieces.push_back(piece);
        h.mNumPoints += n;
    }
    if (h.mPieces.empty())
        throw ReadError("No Piece found");
    return h;
}

/** @brief Stack same-dtype, same-width arrays row-wise; empty when they differ. */
NDArray vtu_concat_rows(std::vector<NDArray>& rParts) {
    if (rParts.size() == 1)
        return std::move(rParts[0]);
    const DType dt = rParts[0].Dtype();
    const std::size_t cols = rParts[0].Shape().size() > 1 ? rParts[0].Shape()[1] : 0;
    std::size_t rows = 0;
    for (const NDArray& rPart : rParts) {
        const std::size_t c = rPart.Shape().size() > 1 ? rPart.Shape()[1] : 0;
        if (rPart.Dtype() != dt || c != cols)
            return NDArray();
        rows += rPart.Shape().empty() ? 0 : rPart.Shape()[0];
    }
    std::vector<std::size_t> shape{rows};
    if (cols)
        shape.push_back(cols);
    NDArray out(dt, shape);
    char* dst = reinterpret_cast<char*>(out.Data());
    for (const NDArray& rPart : rParts) {
        const std::size_t nbytes = rPart.Size() * dtype_size(dt);
        if (nbytes)
            std::memcpy(dst, rPart.Data(), nbytes);
        dst += nbytes;
    }
    return out;
}

/** @brief `<DataArray>` `Name` attributes under @p rSection present in every
 * piece (a multi-piece read drops the others), sorted. */
std::vector<std::string> vtu_array_names(const std::vector<pugi::xml_node>& rPieces,
                                         const char* pSection) {
    std::vector<std::string> names;
    for (std::size_t k = 0; k < rPieces.size(); ++k) {
        std::vector<std::string> here;
        for (pugi::xml_node da : rPieces[k].child(pSection).children("DataArray"))
            here.emplace_back(da.attribute("Name").as_string());
        // The uniform mesh API hands back sorted names; match it so a summary
        // and a real read report data arrays in the same order.
        std::sort(here.begin(), here.end());
        if (k == 0) {
            names = std::move(here);
        } else {
            std::vector<std::string> both;
            std::set_intersection(names.begin(), names.end(), here.begin(), here.end(),
                                  std::back_inserter(both));
            names = std::move(both);
        }
    }
    return names;
}

/** @brief Whether a `<DataArray type=>` is one of the ten numeric types meshio++ holds. */
bool vtu_is_numeric_type(const std::string& rType) {
    return rType == "Float32" || rType == "Float64" || rType == "Int8" || rType == "Int16" ||
           rType == "Int32" || rType == "Int64" || rType == "UInt8" || rType == "UInt16" ||
           rType == "UInt32" || rType == "UInt64";
}

/**
 * @brief Read the `<FieldData>` arrays under @p rNode into `mesh.field_data`.
 *
 * Field data belongs to the dataset, not to a piece: VTK writes it on the
 * `<UnstructuredGrid>` element, before the `<Piece>`, and also accepts it inside one, so the
 * reader looks at both (the piece's overriding the grid's, since `AddFieldData`
 * is insert-or-assign). A non-numeric array (`type="String"`, `"Bit"`) has no
 * meshio++ dtype: it is skipped with a warning rather than failing a read that
 * used to succeed by ignoring the whole section.
 */
void vtu_read_field_data(const pugi::xml_node& rNode, const VtuContext& rCtx,
                         const ReadOptions& rOpts, bool WantData,
                         std::vector<std::pair<std::string, NDArray>>& rOut) {
    for (pugi::xml_node da : rNode.child("FieldData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        // Region arrays are topology, not data: read whatever the options
        // narrow (detail/region_field_data.hpp).
        if (!detail::is_region_field_name(name) && (!WantData || !rOpts.WantsArray(name)))
            continue;
        if (!vtu_is_numeric_type(da.attribute("type").as_string())) {
            log::warn(
                "meshio++: VTU: skipping <FieldData> array '{}' of type '{}' (only numeric "
                "arrays are read)",
                name, da.attribute("type").as_string());
            continue;
        }
        int nc = 0;
        NDArray arr = vtu_read_data_array(da, rCtx, nc);
        if (nc > 1)
            arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
        rOut.emplace_back(name, std::move(arr));
    }
}

/**
 * @brief The field-data names a real read would return: the numeric arrays of the
 * grid's and the piece's `<FieldData>`, sorted and unique.
 *
 * Numeric only, like `vtu_read_field_data`, so a summary never names an array the
 * read skips.
 */
std::vector<std::string> vtu_field_data_names(const vtu_header& rHeader) {
    std::vector<pugi::xml_node> nodes{rHeader.mGrid};
    nodes.insert(nodes.end(), rHeader.mPieces.begin(), rHeader.mPieces.end());
    std::vector<std::string> names;
    for (const pugi::xml_node& rNode : nodes)
        for (pugi::xml_node da : rNode.child("FieldData").children("DataArray"))
            if (vtu_is_numeric_type(da.attribute("type").as_string()) &&
                !detail::is_region_field_name(da.attribute("Name").as_string()))
                names.emplace_back(da.attribute("Name").as_string());
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

}  // namespace

Mesh read_vtu(const std::string& rPath, const ReadOptions& rOpts) {
    VtuSource source;
    vtu_load(rPath, pugi::parse_default, source);
    const vtu_header h = vtu_parse_header(source);
    const VtuContext& ctx = h.mCtx;
    const bool want_data = rOpts.WantsAnyData();
    const bool many = h.mPieces.size() > 1;

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    // VTU's polyhedral stream; empty when the file has none.
    std::vector<std::int64_t> faces, face_offsets;
    std::vector<NDArray> point_parts;
    // Per name, one array per piece; a name some piece lacks is dropped.
    std::map<std::string, std::vector<NDArray>> point_data, cell_data;
    std::int64_t point_base = 0, conn_base = 0, face_base = 0;

    for (std::size_t k = 0; k < h.mPieces.size(); ++k) {
        const pugi::xml_node piece = h.mPieces[k];
        const std::size_t num_points =
            static_cast<std::size_t>(piece.attribute("NumberOfPoints").as_ullong());
        const std::size_t num_cells =
            static_cast<std::size_t>(piece.attribute("NumberOfCells").as_ullong());
        std::vector<std::int64_t> p_conn, p_offsets, p_types, p_faces, p_face_offsets;
        for (pugi::xml_node child : piece.children()) {
            std::string tag = child.name();
            if (tag == "Points") {
                pugi::xml_node da = child.child("DataArray");
                int nc = 0;
                NDArray pts = vtu_read_data_array(da, ctx, nc);
                if (nc <= 0)
                    nc = 3;
                pts.Reshape({num_points, static_cast<std::size_t>(nc)});
                point_parts.push_back(std::move(pts));
            } else if (tag == "Cells") {
                for (pugi::xml_node da : child.children("DataArray")) {
                    int nc = 0;
                    std::string name = da.attribute("Name").as_string();
                    if (name != "connectivity" && name != "offsets" && name != "types" &&
                        name != "faces" && name != "faceoffsets")
                        continue;
                    NDArray arr = vtu_read_data_array(da, ctx, nc);
                    if (name == "connectivity")
                        p_conn = vtu_to_int64(arr);
                    else if (name == "offsets")
                        p_offsets = vtu_to_int64(arr);
                    else if (name == "types")
                        p_types = vtu_to_int64(arr);
                    else if (name == "faces")
                        p_faces = vtu_to_int64(arr);
                    else
                        p_face_offsets = vtu_to_int64(arr);
                }
            } else if (tag == "PointData" || tag == "CellData") {
                if (!want_data)
                    continue;
                auto& rTarget = tag == "PointData" ? point_data : cell_data;
                for (pugi::xml_node da : child.children("DataArray")) {
                    int nc = 0;
                    std::string name = da.attribute("Name").as_string();
                    // The Name attribute is available before the body is
                    // touched, so an unwanted array costs nothing but the
                    // attribute read -- no base64 decode, no inflate, no
                    // allocation.
                    if (!rOpts.WantsArray(name))
                        continue;
                    NDArray arr = vtu_read_data_array(da, ctx, nc);
                    if (nc > 1)
                        arr.Reshape({arr.Size() / nc, static_cast<std::size_t>(nc)});
                    auto& rParts = rTarget[name];
                    if (rParts.size() == k)  // one per piece; a repeat keeps the first
                        rParts.push_back(std::move(arr));
                }
            }
        }
        if (many && (p_offsets.size() != num_cells || p_types.size() != num_cells))
            throw ReadError("VTU: piece " + std::to_string(k) + " has inconsistent cells");
        // One stream across pieces: node ids shift by the points before the
        // piece, offsets by the connectivity before it, face offsets by the
        // face stream before it (-1 marks a cell that is not a polyhedron).
        if (point_base != 0)
            for (auto& rV : p_conn)
                rV += point_base;
        if (conn_base != 0)
            for (auto& rV : p_offsets)
                rV += conn_base;
        for (std::size_t i = 0; i < p_faces.size();) {
            const std::int64_t num_faces = p_faces[i++];
            for (std::int64_t f = 0; f < num_faces && i < p_faces.size(); ++f) {
                const std::int64_t n = p_faces[i++];
                for (std::int64_t j = 0; j < n && i < p_faces.size(); ++j)
                    p_faces[i++] += point_base;
            }
        }
        if (p_face_offsets.empty() && (!faces.empty() || !face_offsets.empty()))
            p_face_offsets.assign(p_types.size(), -1);
        if (!p_face_offsets.empty() && face_offsets.empty() && !types.empty())
            face_offsets.assign(types.size(), -1);
        for (auto& rV : p_face_offsets)
            if (rV >= 0)
                rV += face_base;
        point_base += static_cast<std::int64_t>(num_points);
        conn_base += static_cast<std::int64_t>(p_conn.size());
        face_base += static_cast<std::int64_t>(p_faces.size());
        // The first (usually the only) piece moves in; later ones append.
        const auto append = [](std::vector<std::int64_t>& rAll, std::vector<std::int64_t>& rPart) {
            if (rAll.empty())
                rAll = std::move(rPart);
            else
                rAll.insert(rAll.end(), rPart.begin(), rPart.end());
        };
        append(conn, p_conn);
        append(offsets, p_offsets);
        append(types, p_types);
        append(faces, p_faces);
        append(face_offsets, p_face_offsets);
    }

    if (!point_parts.empty()) {
        if (point_parts.size() != h.mPieces.size())
            throw ReadError("VTU: a piece has no Points");
        NDArray points = vtu_concat_rows(point_parts);
        if (points.Size() == 0 && h.mNumPoints > 0)
            throw ReadError("VTU: pieces disagree on the point dtype or dimension");
        mesh.AssignPoints(std::move(points));
    }

    auto merged = [&](std::map<std::string, std::vector<NDArray>>& rData, const char* pWhat,
                      auto&& rAdd) {
        for (auto& [name, parts] : rData) {
            NDArray arr = parts.size() == h.mPieces.size() ? vtu_concat_rows(parts) : NDArray();
            if (parts.size() != h.mPieces.size() ||
                (arr.Size() == 0 && !parts.empty() && parts[0].Size() != 0)) {
                log::warn(
                    "VTU: {} data '{}' is missing from, or differs between, pieces; "
                    "dropped",
                    pWhat, name);
                continue;
            }
            rAdd(name, std::move(arr));
        }
    };
    std::unordered_map<std::string, NDArray> cell_data_raw;
    merged(point_data, "point", [&](const std::string& rName, NDArray&& rArr) {
        mesh.AddPointData(rName, std::move(rArr));
    });
    merged(cell_data, "cell", [&](const std::string& rName, NDArray&& rArr) {
        cell_data_raw.emplace(rName, std::move(rArr));
    });

    std::vector<std::pair<std::string, NDArray>> field_arrays;
    vtu_read_field_data(h.mGrid, ctx, rOpts, want_data, field_arrays);
    for (const pugi::xml_node& rPiece : h.mPieces)
        vtu_read_field_data(rPiece, ctx, rOpts, want_data, field_arrays);

    detail::check_vtk_cell_arrays(conn.size(), offsets, types, cell_data_raw);
    std::vector<std::int64_t> file_to_global;
    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw,
                              faces.empty() ? nullptr : &faces, face_offsets, mesh,
                              &file_to_global);
    detail::regions_from_field_arrays(mesh, field_arrays, &file_to_global, "vtu");
    return mesh;
}

MeshMetadata read_vtu_metadata(const std::string& rPath, const ReadOptions&) {
    // parse_minimal skips escape expansion and EOL normalization over the
    // base64 bodies. It does NOT avoid reading the file: pugixml always
    // materializes PCDATA. The real saving below is skipping base64 decode,
    // decompression, allocation and byte-swapping for every array we don't
    // touch -- a solid multiple, not an asymptotic change. See read_options.hpp.
    VtuSource source;
    vtu_load(rPath, pugi::parse_minimal, source);
    const vtu_header h = vtu_parse_header(source);

    MeshMetadata meta;
    meta.mNumPoints = h.mNumPoints;  // attributes -- free

    // Point dimension comes from the Points DataArray's NumberOfComponents
    // attribute, so the coordinates themselves are never decoded.
    pugi::xml_node points_da = h.mPieces[0].child("Points").child("DataArray");
    const int point_nc = points_da ? points_da.attribute("NumberOfComponents").as_int(0) : 0;
    meta.mPointDim = point_nc > 0 ? static_cast<std::size_t>(point_nc) : 3;

    // 'types' is one value per cell and is the only array a summary must decode;
    // 'offsets' is additionally needed only when a variable-node-count type
    // (polygon, VTK_LAGRANGE_*) is present.
    auto piece_array = [&](const pugi::xml_node& rPiece, const char* pName) {
        for (pugi::xml_node da : rPiece.child("Cells").children("DataArray"))
            if (std::string(da.attribute("Name").as_string()) == pName) {
                int nc = 0;
                return vtu_to_int64(vtu_read_data_array(da, h.mCtx, nc));
            }
        return std::vector<std::int64_t>();
    };
    std::vector<std::int64_t> types, offsets;
    for (const pugi::xml_node& rPiece : h.mPieces) {
        std::vector<std::int64_t> t = piece_array(rPiece, "types");
        types.insert(types.end(), t.begin(), t.end());
    }
    if (detail::cells_need_offsets(types)) {
        std::int64_t base = 0;
        for (const pugi::xml_node& rPiece : h.mPieces) {
            std::vector<std::int64_t> o = piece_array(rPiece, "offsets");
            for (auto& rV : o)
                rV += base;
            if (!o.empty())
                base = o.back();
            offsets.insert(offsets.end(), o.begin(), o.end());
        }
    }
    meta.mCellBlocks = detail::summarize_cells(offsets, types);

    meta.mPointDataNames = vtu_array_names(h.mPieces, "PointData");
    meta.mCellDataNames = vtu_array_names(h.mPieces, "CellData");
    meta.mFieldDataNames = vtu_field_data_names(h);
    {
        std::vector<std::pair<std::string, std::size_t>> arrays;
        for (const pugi::xml_node& rNode : {h.mGrid, h.mPieces[0]})
            for (pugi::xml_node da : rNode.child("FieldData").children("DataArray"))
                arrays.emplace_back(da.attribute("Name").as_string(),
                                    da.attribute("NumberOfTuples").as_ullong(0));
        meta.mRegions = detail::region_summaries_from_field_names(arrays);
    }

    // No bounding box: it would require decoding the point coordinates, which
    // are usually the largest array in the file -- exactly what this path exists
    // to avoid. metadata_from_mesh fills it in on the full-read fallback.
    meta.mHasBBox = false;
    return meta;
}

}  // namespace meshioplusplus
