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
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/frd.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/sym3_eigen.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

namespace {

constexpr std::size_t frd_values_per_line = 6;
constexpr std::size_t frd_value_width = 12;

// FRD element type -> cell type, node count and permutation. `connectivity[k] =
// frd_nodes[permutation[k]]`; a null permutation is the identity. The he20 and pe15
// mid-edge groups and the be3 mid node sit elsewhere than in Abaqus order (confirmed
// against ccx 2.23 output for the same `.inp`); types 7-10 are cgx's own shells, which
// ccx expands into solids and never writes.
struct FrdTypeSpec {
    const char* mName;
    std::size_t mNodes;
    const int* mPermutation;
};

constexpr int frd_perm_he20[20] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                   10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int frd_perm_pe15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};
constexpr int frd_perm_be3[3] = {0, 2, 1};

const FrdTypeSpec* frd_type_spec(int Type) {
    static const std::array<FrdTypeSpec, 12> table = {{
        {"hexahedron", 8, nullptr},
        {"wedge", 6, nullptr},
        {"tetra", 4, nullptr},
        {"hexahedron20", 20, frd_perm_he20},
        {"wedge15", 15, frd_perm_pe15},
        {"tetra10", 10, nullptr},
        {"triangle", 3, nullptr},
        {"triangle6", 6, nullptr},
        {"quad", 4, nullptr},
        {"quad8", 8, nullptr},
        {"line", 2, nullptr},
        {"line3", 3, frd_perm_be3},
    }};
    if (Type < 1 || Type > 12)
        return nullptr;
    return &table[static_cast<std::size_t>(Type - 1)];
}

bool frd_is_tensor_name(const std::string& rName) {
    return rName == "STRESS" || rName == "TOSTRAIN" || rName == "MESTRAIN" || rName == "ZZSTR";
}

[[noreturn]] void frd_fail(const std::string& rMessage) {
    throw ReadError("CalculiX FRD: " + rMessage);
}

bool frd_starts_with(std::string_view Line, std::string_view Prefix) {
    return Line.size() >= Prefix.size() && Line.compare(0, Prefix.size(), Prefix) == 0;
}

/// Columns [Pos, Pos + Len) of @p Line, clipped to its end.
std::string_view frd_field(std::string_view Line, std::size_t Pos, std::size_t Len) {
    if (Pos >= Line.size())
        return {};
    return Line.substr(Pos, Len);
}

std::string_view frd_strip(std::string_view Text) {
    std::size_t b = 0;
    std::size_t e = Text.size();
    while (b < e && (Text[b] == ' ' || Text[b] == '\t'))
        ++b;
    while (e > b && (Text[e - 1] == ' ' || Text[e - 1] == '\t' || Text[e - 1] == '\r'))
        --e;
    return Text.substr(b, e - b);
}

std::int64_t frd_int(std::string_view Text, const char* pWhere) {
    const std::string_view s = frd_strip(Text);
    if (s.empty())
        return 0;
    std::size_t i = 0;
    bool negative = false;
    if (s[0] == '-' || s[0] == '+') {
        negative = s[0] == '-';
        i = 1;
    }
    if (i >= s.size())
        frd_fail("invalid integer field '" + std::string(s) + "' in " + pWhere);
    std::int64_t value = 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9')
            frd_fail("invalid integer field '" + std::string(s) + "' in " + pWhere);
        value = value * 10 + (s[i] - '0');
    }
    return negative ? -value : value;
}

double frd_real(std::string_view Text, const char* pWhere) {
    const std::string_view s = frd_strip(Text);
    char buf[64];
    if (s.empty() || s.size() >= sizeof(buf))
        frd_fail("invalid real field '" + std::string(s) + "' in " + pWhere);
    for (std::size_t i = 0; i < s.size(); ++i)
        buf[i] = s[i] == 'D' ? 'E' : (s[i] == 'd' ? 'e' : s[i]);
    buf[s.size()] = '\0';
    const char* end = nullptr;
    const double value = detail::parse_double(buf, end);
    if (end != buf + s.size())
        frd_fail("invalid real field '" + std::string(s) + "' in " + pWhere);
    return value;
}

/// Trailing format flag of a `2C`/`3C` header line (0 short, 1 long).
int frd_flag(std::string_view Line, int Default) {
    std::string_view last;
    std::size_t tokens = 0;
    std::size_t i = 0;
    while (i < Line.size()) {
        while (i < Line.size() && (Line[i] == ' ' || Line[i] == '\t' || Line[i] == '\r'))
            ++i;
        const std::size_t start = i;
        while (i < Line.size() && Line[i] != ' ' && Line[i] != '\t' && Line[i] != '\r')
            ++i;
        if (i > start) {
            last = Line.substr(start, i - start);
            ++tokens;
        }
    }
    if (tokens >= 3) {
        bool digits = !last.empty();
        for (char c : last)
            digits = digits && c >= '0' && c <= '9';
        if (digits)
            return static_cast<int>(frd_int(last, "a header record"));
    }
    return Default;
}

std::size_t frd_width(int Flag) {
    if (Flag == 2)
        frd_fail("binary .frd files are not supported");
    return Flag == 0 ? 5 : 10;
}

struct FrdBlock {
    std::string mName;
    std::size_t mDataComps = 0;
    std::size_t mWidth = 10;
    std::size_t mFirst = 0;
    std::size_t mLast = 0;
};

struct FrdFrame {
    std::string mKey;
    std::string mValueText;
    double mValue = 0.0;
    std::int64_t mAnalysis = 0;
    std::int64_t mStep = 0;
    std::vector<FrdBlock> mBlocks;
};

struct FrdElement {
    int mType = 0;
    std::int64_t mGroup = 0;
    std::int64_t mMaterial = 0;
    std::size_t mFirstNode = 0;  // into FrdFile::mElementNodes
    std::size_t mNumNodes = 0;
};

class FrdFile {
public:
    explicit FrdFile(const std::string& rPath) {
        auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
        if (!in)
            frd_fail("could not read " + rPath);
        in.seekg(0, std::ios::end);
        const std::streamoff size = in.tellg();
        in.seekg(0, std::ios::beg);
        mText.resize(size > 0 ? static_cast<std::size_t>(size) : 0);
        if (!mText.empty())
            in.read(mText.data(), static_cast<std::streamsize>(mText.size()));
        IndexLines();
        Parse();
    }

    std::vector<std::int64_t> mNodeIds;
    std::vector<double> mCoords;
    std::vector<FrdElement> mElements;
    std::vector<std::int64_t> mElementNodes;
    std::vector<FrdFrame> mFrames;

    std::unordered_map<std::int64_t, std::int64_t> NodeIndex() const {
        std::unordered_map<std::int64_t, std::int64_t> index;
        index.reserve(mNodeIds.size() * 2);
        for (std::size_t i = 0; i < mNodeIds.size(); ++i)
            index.emplace(mNodeIds[i], static_cast<std::int64_t>(i));
        return index;
    }

    /// Parse @p rBlock into a (num nodes x DataComps) row-major array, NaN where absent.
    void ReadBlock(const FrdBlock& rBlock,
                   const std::unordered_map<std::int64_t, std::int64_t>& rIndex,
                   std::vector<double>& rOut) const {
        const std::size_t nc = rBlock.mDataComps;
        rOut.assign(mNodeIds.size() * nc, std::numeric_limits<double>::quiet_NaN());
        const std::size_t end = 3 + rBlock.mWidth;
        std::size_t pos = rBlock.mFirst;
        while (pos < rBlock.mLast) {
            std::string_view line = mLines[pos++];
            if (!frd_starts_with(line, " -1"))
                continue;
            const std::int64_t node = frd_int(frd_field(line, 3, rBlock.mWidth), "a result line");
            const auto it = rIndex.find(node);
            if (it == rIndex.end())
                frd_fail("result for " + rBlock.mName + " refers to undefined node " +
                         std::to_string(node));
            double* row = rOut.data() + static_cast<std::size_t>(it->second) * nc;
            std::size_t got = 0;
            while (true) {
                const std::size_t take = std::min(frd_values_per_line, nc - got);
                for (std::size_t k = 0; k < take; ++k) {
                    const std::string_view field =
                        frd_field(line, end + k * frd_value_width, frd_value_width);
                    if (frd_strip(field).empty())
                        frd_fail("result line for " + rBlock.mName + " is missing values");
                    row[got + k] = frd_real(field, "a result line");
                }
                got += take;
                if (got >= nc)
                    break;
                if (pos >= rBlock.mLast || !frd_starts_with(mLines[pos], " -2"))
                    frd_fail("result line for " + rBlock.mName + " is missing values");
                line = mLines[pos++];
            }
        }
    }

private:
    std::string mText;
    std::vector<std::string_view> mLines;

    void IndexLines() {
        const std::string_view all(mText);
        std::size_t start = 0;
        while (start <= all.size()) {
            std::size_t stop = all.find('\n', start);
            if (stop == std::string_view::npos)
                stop = all.size();
            std::string_view line = all.substr(start, stop - start);
            while (!line.empty() && line.back() == '\r')
                line.remove_suffix(1);
            mLines.push_back(line);
            start = stop + 1;
        }
    }

    void Parse() {
        const std::size_t n = mLines.size();
        std::size_t pos = 0;
        int flag = 1;
        bool seen_nodes = false;
        bool seen_elements = false;
        while (pos < n) {
            const std::string_view line = mLines[pos++];
            if (frd_starts_with(line, "    2C")) {
                flag = frd_flag(line, flag);
                if (seen_nodes) {
                    log::warn("{}", "CalculiX FRD: a second node block was ignored");
                    pos = SkipBlock(pos);
                } else {
                    pos = ReadNodes(pos, frd_width(flag));
                    seen_nodes = true;
                }
            } else if (frd_starts_with(line, "    3C")) {
                flag = frd_flag(line, flag);
                if (seen_elements) {
                    log::warn("{}", "CalculiX FRD: a second element block was ignored");
                    pos = SkipBlock(pos);
                } else {
                    pos = ReadElements(pos, frd_width(flag));
                    seen_elements = true;
                }
            } else if (frd_starts_with(line, "  100C")) {
                pos = ReadFrame(pos, line, flag);
            } else if (frd_starts_with(line, "  9999")) {
                break;
            }
        }
        if (!seen_nodes)
            frd_fail("no node block (2C record): not a CalculiX result file");
    }

    std::size_t SkipBlock(std::size_t Pos) const {
        while (Pos < mLines.size() && !frd_starts_with(mLines[Pos], " -3"))
            ++Pos;
        return Pos + 1;
    }

    std::size_t ReadNodes(std::size_t Pos, std::size_t W) {
        const std::size_t end = 3 + W;
        while (Pos < mLines.size()) {
            const std::string_view line = mLines[Pos++];
            if (frd_starts_with(line, " -3"))
                break;
            if (!frd_starts_with(line, " -1"))
                continue;
            mNodeIds.push_back(frd_int(frd_field(line, 3, W), "a node line"));
            for (std::size_t k = 0; k < 3; ++k)
                mCoords.push_back(frd_real(
                    frd_field(line, end + k * frd_value_width, frd_value_width), "a node line"));
        }
        return Pos;
    }

    std::size_t ReadElements(std::size_t Pos, std::size_t W) {
        const std::size_t end = 3 + W;
        bool have_current = false;
        while (Pos < mLines.size()) {
            const std::string_view line = mLines[Pos++];
            if (frd_starts_with(line, " -3"))
                break;
            if (frd_starts_with(line, " -1")) {
                FrdElement el;
                el.mType = static_cast<int>(frd_int(frd_field(line, end, 5), "an element line"));
                el.mGroup = frd_int(frd_field(line, end + 5, 5), "an element line");
                el.mMaterial = frd_int(frd_field(line, end + 10, 5), "an element line");
                el.mFirstNode = mElementNodes.size();
                mElements.push_back(el);
                have_current = true;
            } else if (frd_starts_with(line, " -2") && have_current) {
                std::string_view body = line.substr(3);
                while (!body.empty() &&
                       (body.back() == ' ' || body.back() == '\t' || body.back() == '\r'))
                    body.remove_suffix(1);
                if (body.size() % W)
                    frd_fail("element line '" + std::string(line) +
                             "' is not a whole number of ids");
                for (std::size_t k = 0; k < body.size(); k += W) {
                    mElementNodes.push_back(frd_int(body.substr(k, W), "an element line"));
                    ++mElements.back().mNumNodes;
                }
            }
        }
        return Pos;
    }

    std::size_t ReadFrame(std::size_t Pos, std::string_view Header, int Flag) {
        int frame_flag = Flag;
        if (Header.size() >= 75)
            frame_flag = static_cast<int>(frd_int(frd_field(Header, 73, 2), "a 100C record"));
        const std::size_t w = frd_width(frame_flag);
        const std::string key(frd_field(Header, 6, 6));
        const std::string value_text(frd_field(Header, 12, 12));
        FrdFrame* frame = nullptr;
        for (FrdFrame& f : mFrames)
            if (f.mKey == key && f.mValueText == value_text) {
                frame = &f;
                break;
            }
        if (frame == nullptr) {
            FrdFrame f;
            f.mKey = key;
            f.mValueText = value_text;
            f.mValue = frd_real(value_text, "a 100C record");
            f.mAnalysis = frd_int(frd_field(Header, 56, 2), "a 100C record");
            f.mStep = frd_int(frd_field(Header, 58, 5), "a 100C record");
            mFrames.push_back(std::move(f));
            frame = &mFrames.back();
        }
        while (Pos < mLines.size()) {
            const std::string_view line = mLines[Pos++];
            if (frd_starts_with(line, " -4")) {
                FrdBlock block;
                block.mName = std::string(frd_strip(frd_field(line, 5, 8)));
                const std::int64_t ncomps = frd_int(frd_field(line, 13, 5), "a -4 record");
                block.mWidth = w;
                std::int64_t calculated = 0;
                while (Pos < mLines.size() && frd_starts_with(mLines[Pos], " -5")) {
                    if (frd_int(frd_field(mLines[Pos], 33, 5), "a -5 record") == 1)
                        ++calculated;
                    ++Pos;
                }
                block.mDataComps =
                    static_cast<std::size_t>(std::max<std::int64_t>(ncomps - calculated, 0));
                block.mFirst = Pos;
                while (Pos < mLines.size() && !frd_starts_with(mLines[Pos], " -3"))
                    ++Pos;
                block.mLast = Pos;
                ++Pos;
                frame->mBlocks.push_back(std::move(block));
                break;
            }
            if (frd_starts_with(line, "  100") || frd_starts_with(line, "    1P") ||
                frd_starts_with(line, "  9999")) {
                --Pos;
                break;
            }
        }
        return Pos;
    }
};

NDArray frd_scalar_array(DType Type, double Value) {
    NDArray out(Type, {std::size_t{1}});
    if (Type == DType::Float64)
        *out.As<double>() = Value;
    else
        *out.As<std::int64_t>() = static_cast<std::int64_t>(Value);
    return out;
}

}  // namespace

Mesh read_frd(const std::string& rPath, const ReadOptions& rOpts) {
    return read_frd(rPath, rOpts, FrdReadOptions{});
}

Mesh read_frd(const std::string& rPath, const ReadOptions& rOpts, const FrdReadOptions& rFrdOpts) {
    const FrdFile file(rPath);
    const auto index = file.NodeIndex();

    Mesh mesh;
    const std::size_t npts = file.mNodeIds.size();
    NDArray points(DType::Float64, {npts, std::size_t{3}});
    std::copy(file.mCoords.begin(), file.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));

    // One block per contiguous run of one cell type.
    struct Run {
        const FrdTypeSpec* mSpec;
        std::vector<std::int64_t> mConn;
        std::vector<std::int64_t> mGroup;
        std::vector<std::int64_t> mMaterial;
    };
    std::vector<Run> runs;
    std::vector<int> skipped;
    for (const FrdElement& el : file.mElements) {
        const FrdTypeSpec* spec = frd_type_spec(el.mType);
        if (spec == nullptr) {
            if (std::find(skipped.begin(), skipped.end(), el.mType) == skipped.end())
                skipped.push_back(el.mType);
            continue;
        }
        if (el.mNumNodes != spec->mNodes)
            frd_fail("element type " + std::to_string(el.mType) + " (" + spec->mName + ") needs " +
                     std::to_string(spec->mNodes) + " nodes, found " +
                     std::to_string(el.mNumNodes));
        if (runs.empty() || runs.back().mSpec != spec)
            runs.push_back(Run{spec, {}, {}, {}});
        Run& run = runs.back();
        const std::size_t base = run.mConn.size();
        run.mConn.resize(base + spec->mNodes);
        for (std::size_t k = 0; k < spec->mNodes; ++k) {
            const std::size_t src =
                spec->mPermutation ? static_cast<std::size_t>(spec->mPermutation[k]) : k;
            const std::int64_t id = file.mElementNodes[el.mFirstNode + src];
            const auto it = index.find(id);
            if (it == index.end())
                frd_fail("an element references undefined node " + std::to_string(id));
            run.mConn[base + k] = it->second;
        }
        run.mGroup.push_back(el.mGroup);
        run.mMaterial.push_back(el.mMaterial);
    }
    if (!skipped.empty()) {
        std::sort(skipped.begin(), skipped.end());
        std::string list;
        for (int t : skipped)
            list += (list.empty() ? "" : ", ") + std::to_string(t);
        log::warn("{}", "CalculiX FRD: skipped elements of unsupported type(s) " + list);
    }
    std::vector<NDArray> group_blocks;
    std::vector<NDArray> material_blocks;
    for (Run& run : runs) {
        const std::size_t count = run.mGroup.size();
        NDArray conn(DType::Int64, {count, run.mSpec->mNodes});
        std::copy(run.mConn.begin(), run.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(run.mSpec->mName, std::move(conn));
        NDArray group(DType::Int64, {count});
        std::copy(run.mGroup.begin(), run.mGroup.end(), group.As<std::int64_t>());
        group_blocks.push_back(std::move(group));
        NDArray material(DType::Int64, {count});
        std::copy(run.mMaterial.begin(), run.mMaterial.end(), material.As<std::int64_t>());
        material_blocks.push_back(std::move(material));
    }
    if (!runs.empty()) {
        mesh.AddCellData("frd:group", std::move(group_blocks));
        mesh.AddCellData("frd:material", std::move(material_blocks));
    }

    if (file.mFrames.empty()) {
        rOpts.ResolveTimeStep(0);  // refuses anything but the first/last "step" of no steps
        return mesh;
    }
    const FrdFrame& frame = file.mFrames[rOpts.ResolveTimeStep(file.mFrames.size())];
    mesh.AddFieldData(kSequenceTimeKey, frd_scalar_array(DType::Float64, frame.mValue));
    mesh.AddFieldData("frd:step", frd_scalar_array(DType::Int64, static_cast<double>(frame.mStep)));
    mesh.AddFieldData("frd:analysis",
                      frd_scalar_array(DType::Int64, static_cast<double>(frame.mAnalysis)));
    if (!rOpts.WantsAnyData())
        return mesh;

    std::vector<std::string> used;
    std::vector<double> values;
    for (const FrdBlock& block : frame.mBlocks) {
        std::string name = block.mName;
        for (int k = 2; std::find(used.begin(), used.end(), name) != used.end(); ++k)
            name = block.mName + "_" + std::to_string(k);
        used.push_back(name);
        const bool tensor =
            rFrdOpts.mDerived && frd_is_tensor_name(block.mName) && block.mDataComps == 6;
        const bool want_mises = tensor && rOpts.WantsArray(name + "_mises");
        const bool want_principal = tensor && rOpts.WantsArray(name + "_principal");
        const bool want_raw = rOpts.WantsArray(name);
        if (!(want_raw || want_mises || want_principal))
            continue;
        file.ReadBlock(block, index, values);
        const std::size_t nc = block.mDataComps;
        if (want_raw) {
            NDArray data =
                nc == 1 ? NDArray(DType::Float64, {npts}) : NDArray(DType::Float64, {npts, nc});
            std::copy(values.begin(), values.end(), data.As<double>());
            mesh.AddPointData(name, std::move(data));
        }
        if (want_mises) {
            NDArray data(DType::Float64, {npts});
            double* out = data.As<double>();
            for (std::size_t i = 0; i < npts; ++i)
                out[i] = detail::sym3_mises(values.data() + i * 6);
            mesh.AddPointData(name + "_mises", std::move(data));
        }
        if (want_principal) {
            NDArray data(DType::Float64, {npts, std::size_t{3}});
            double* out = data.As<double>();
            for (std::size_t i = 0; i < npts; ++i) {
                const double* t = values.data() + i * 6;
                bool finite = true;
                for (int k = 0; k < 6; ++k)
                    finite = finite && std::isfinite(t[k]);
                if (finite)
                    detail::sym3_principal(t, out + i * 3);
                else
                    std::fill(out + i * 3, out + i * 3 + 3,
                              std::numeric_limits<double>::quiet_NaN());
            }
            mesh.AddPointData(name + "_principal", std::move(data));
        }
    }
    return mesh;
}

MeshMetadata read_frd_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    const FrdFile file(rPath);
    // No header-only path: the same full-read-plus-override shape EnSight's has.
    MeshMetadata meta = metadata_from_mesh(read_frd(rPath, ReadOptions{}));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "frd";
    meta.mTimeValues.clear();
    for (const FrdFrame& frame : file.mFrames)
        meta.mTimeValues.push_back(frame.mValue);
    return meta;
}

}  // namespace meshioplusplus
