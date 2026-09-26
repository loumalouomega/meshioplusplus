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
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "vtk_preflight.hpp"
#include "../detail/vtu_decode.hpp"

namespace meshioplusplus {

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
    for (const std::uint64_t z : sizes)
        comp_total += std::min<std::uint64_t>(z, std::uint64_t{1} << 40);
    const std::uint64_t ceiling = (std::uint64_t{1} << 16) * (comp_total + 1) + (1u << 20);
    if (num_blocks > 0 && max_block <= ceiling && last_block <= ceiling &&
        num_blocks - 1 <= ceiling / std::max<std::uint64_t>(max_block, 1)) {
        const std::uint64_t want = (num_blocks - 1) * max_block + last_block;
        if (want <= ceiling)
            out.reserve(static_cast<std::size_t>(want));
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
        out.insert(out.end(), dec.begin(), dec.begin() + std::min(dec.size(), expected));
    }
    return out;
}

NDArray vtu_array_from_bytes(const std::vector<unsigned char>& rBytes, DType Dt, bool BigEndian) {
    const std::size_t isz = dtype_size(Dt);
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

/** @brief How the arrays of one file are framed: codec, header width, byte
 * order and, when present, the `<AppendedData>` payload. */
struct VtuContext {
    detail::VtkCodec mCodec = detail::VtkCodec::None;
    std::size_t mHeaderSize = 4;
    bool mBigEndian = false;
    // Raw appended payload (bytes after the opening '_'), or base64 text.
    const unsigned char* mRaw = nullptr;
    std::size_t mRawLen = 0;
    const char* mBase64 = nullptr;
    std::size_t mBase64Len = 0;
};

NDArray vtu_read_data_array(const pugi::xml_node& rDa, const VtuContext& rCtx,
                            int& rNumComponents) {
    std::string fmt = rDa.attribute("format").as_string("ascii");
    DType dt = detail::dtype_from_vtu(rDa.attribute("type").as_string());
    rNumComponents = rDa.attribute("NumberOfComponents").as_int(0);

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
        return vtu_array_from_bytes(vtu_decode_sequential(src, rCtx.mCodec, rCtx.mHeaderSize, true),
                                    dt, true);
    }
    if (fmt == "appended") {
        // strtoull skips the padding some writers put around the offset.
        const std::uint64_t offset =
            std::strtoull(rDa.attribute("offset").as_string("0"), nullptr, 10);
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
}

/**
 * @brief The XML document and the file's bytes when they are needed.
 *
 * A raw `<AppendedData>` payload is not XML text (it may hold any byte, `<`
 * included): the file is read as bytes, the payload cut out, and only the text
 * around it parsed. Every other file parses as before.
 */
struct VtuSource {
    pugi::xml_document mDoc;
    std::string mBytes;
    std::size_t mRawStart = 0;
    std::size_t mRawStop = 0;
    bool mIsRaw = false;
};

void vtu_load(const std::string& rPath, unsigned int ParseOptions, VtuSource& rSource) {
    detail::vtk_preflight(rPath, "UnstructuredGrid",
                          "lzma-compressed VTU not supported by the C++ reader");
    pugi::xml_parse_result res = rSource.mDoc.load_file(rPath.c_str(), ParseOptions);
    if (res) {
        pugi::xml_node app = rSource.mDoc.child("VTKFile").child("AppendedData");
        if (!app || std::string(app.attribute("encoding").as_string("base64")) != "raw")
            return;
        // A raw payload that happened to parse as text still has to be read as bytes.
    }
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    rSource.mBytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    const std::string& b = rSource.mBytes;
    const std::size_t tag = b.find("<AppendedData");
    const std::size_t tag_end = tag == std::string::npos ? tag : b.find('>', tag);
    const bool raw = tag_end != std::string::npos &&
                     b.substr(tag, tag_end - tag).find("\"raw\"") != std::string::npos;
    if (!raw) {
        if (!res)
            throw ReadError(std::string("VTU XML parse failed: ") + res.description());
        return;
    }
    const std::size_t underscore = b.find('_', tag_end);
    const std::size_t stop = b.rfind("</AppendedData>");
    if (underscore == std::string::npos || stop == std::string::npos || stop <= underscore)
        throw ReadError("VTU: AppendedData is not closed");
    const std::string xml = b.substr(0, tag_end + 1) + b.substr(stop);
    rSource.mDoc.reset();
    res = rSource.mDoc.load_buffer(xml.data(), xml.size(), ParseOptions);
    if (!res)
        throw ReadError(std::string("VTU XML parse failed: ") + res.description());
    rSource.mIsRaw = true;
    rSource.mRawStart = underscore + 1;
    rSource.mRawStop = stop;
}

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
    const std::string compressor = root.attribute("compressor").as_string("");
    if (compressor.empty())
        h.mCtx.mCodec = detail::VtkCodec::None;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::Zlib))
        h.mCtx.mCodec = detail::VtkCodec::Zlib;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::LZ4))
        h.mCtx.mCodec = detail::VtkCodec::LZ4;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::ZSTD))
        h.mCtx.mCodec = detail::VtkCodec::ZSTD;
    else if (compressor == detail::vtk_codec_compressor(detail::VtkCodec::LZMA))
        throw ReadError("lzma-compressed VTU not supported by the C++ reader");
    else
        throw ReadError("Unknown VTU compressor '" + compressor + "'");
    // Fail early and actionably when the file needs a codec this build lacks,
    // rather than at the first array body.
    detail::vtk_codec_require_read(h.mCtx.mCodec);

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    h.mCtx.mHeaderSize = (header_type == "UInt64") ? 8 : 4;
    const std::string byte_order = root.attribute("byte_order").as_string("LittleEndian");
    if (byte_order != "LittleEndian" && byte_order != "BigEndian")
        throw ReadError("Unknown VTU byte order '" + byte_order + "'");
    h.mCtx.mBigEndian = byte_order == "BigEndian";

    pugi::xml_node grid = root.child("UnstructuredGrid");
    if (!grid)
        throw ReadError("No UnstructuredGrid found");

    if (rSource.mIsRaw) {
        h.mCtx.mRaw =
            reinterpret_cast<const unsigned char*>(rSource.mBytes.data()) + rSource.mRawStart;
        h.mCtx.mRawLen = rSource.mRawStop - rSource.mRawStart;
    } else if (pugi::xml_node app = root.child("AppendedData")) {
        const std::string encoding = app.attribute("encoding").as_string("base64");
        if (encoding != "base64")
            throw ReadError("Unknown VTU AppendedData encoding '" + encoding + "'");
        const char* text = app.text().get();
        while (*text && std::isspace(static_cast<unsigned char>(*text)))
            ++text;
        if (*text != '_')
            throw ReadError("VTU: AppendedData does not start with '_'");
        h.mCtx.mBase64 = text + 1;
        h.mCtx.mBase64Len = std::strlen(text + 1);
    }

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
                         const ReadOptions& rOpts, Mesh& rMesh) {
    for (pugi::xml_node da : rNode.child("FieldData").children("DataArray")) {
        const std::string name = da.attribute("Name").as_string();
        if (!rOpts.WantsArray(name))
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
        rMesh.AddFieldData(name, std::move(arr));
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
            if (vtu_is_numeric_type(da.attribute("type").as_string()))
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

    if (want_data) {
        vtu_read_field_data(h.mGrid, ctx, rOpts, mesh);
        for (const pugi::xml_node& rPiece : h.mPieces)
            vtu_read_field_data(rPiece, ctx, rOpts, mesh);
    }

    detail::check_vtk_cell_arrays(conn.size(), offsets, types, cell_data_raw);
    detail::reconstruct_cells(conn.data(), offsets, types, cell_data_raw,
                              faces.empty() ? nullptr : &faces, face_offsets, mesh);
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

    // No bounding box: it would require decoding the point coordinates, which
    // are usually the largest array in the file -- exactly what this path exists
    // to avoid. metadata_from_mesh fills it in on the full-read fallback.
    meta.mHasBBox = false;
    return meta;
}

}  // namespace meshioplusplus
