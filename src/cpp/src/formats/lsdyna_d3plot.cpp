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
#include <bit>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <ios>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string_view>
#include <tuple>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/lsdyna_d3plot.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "../detail/typed_view.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

constexpr std::int64_t kD3FemzipNmmat = 76893465;
constexpr double kD3EofMarker = -999999.0;
const double kD3Nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t kD3Npos = std::numeric_limits<std::size_t>::max();

[[noreturn]] void d3_fail(const std::string& rMessage) {
    throw ReadError("LS-DYNA d3plot: " + rMessage);
}

std::int64_t d3_zero_based(std::int64_t Id) {
    // Parts are converted back to user ids later, so reject INT64_MIN
    // rather than letting its wrapped index overflow that inverse conversion.
    if (Id == std::numeric_limits<std::int64_t>::min())
        d3_fail("a one-based index is out of range");
    return detail::zero_based(Id);
}

std::int64_t d3_digit(std::int64_t Value, int I) {
    // Unsigned magnitude: negating INT64_MIN is undefined.
    std::uint64_t v =
        Value < 0 ? 0 - static_cast<std::uint64_t>(Value) : static_cast<std::uint64_t>(Value);
    for (int k = 0; k < I; ++k)
        v /= 10;
    return static_cast<std::int64_t>(v % 10);
}

// --- words ---------------------------------------------------------------------------

// Word access into a byte buffer: 4- or 8-byte words in a given byte order.
struct D3Words {
    const char* mData = nullptr;
    std::size_t mSize = 0;
    int mWs = 4;
    bool mSwap = false;

    std::size_t NumWords() const { return mSize / static_cast<std::size_t>(mWs); }

    void Need(std::size_t Pos, std::size_t N) const {
        // Compared in words, without multiplying: (Pos + N) * mWs can wrap.
        if (Pos > NumWords() || N > NumWords() - Pos)
            d3_fail("the file is truncated");
    }

    std::int64_t IntAt(std::size_t Pos) const {
        const char* p = mData + Pos * static_cast<std::size_t>(mWs);
        if (mWs == 4) {
            std::uint32_t u;
            std::memcpy(&u, p, 4);
            if (mSwap)
                u = detail::bswap32(u);
            return static_cast<std::int32_t>(u);
        }
        std::uint64_t u;
        std::memcpy(&u, p, 8);
        if (mSwap)
            u = detail::bswap64(u);
        return static_cast<std::int64_t>(u);
    }

    double FloatAt(std::size_t Pos) const {
        const char* p = mData + Pos * static_cast<std::size_t>(mWs);
        if (mWs == 4) {
            std::uint32_t u;
            std::memcpy(&u, p, 4);
            if (mSwap)
                u = detail::bswap32(u);
            return static_cast<double>(std::bit_cast<float>(u));
        }
        std::uint64_t u;
        std::memcpy(&u, p, 8);
        if (mSwap)
            u = detail::bswap64(u);
        return std::bit_cast<double>(u);
    }

    std::int64_t Int(std::size_t Pos) const {
        Need(Pos, 1);
        return IntAt(Pos);
    }

    double Float(std::size_t Pos) const {
        Need(Pos, 1);
        return FloatAt(Pos);
    }

    std::vector<std::int64_t> Ints(std::size_t Pos, std::size_t N) const {
        Need(Pos, N);
        std::vector<std::int64_t> out(N);
        for (std::size_t i = 0; i < N; ++i)
            out[i] = IntAt(Pos + i);
        return out;
    }

    std::vector<double> Floats(std::size_t Pos, std::size_t N) const {
        Need(Pos, N);
        std::vector<double> out(N);
        for (std::size_t i = 0; i < N; ++i)
            out[i] = FloatAt(Pos + i);
        return out;
    }

    // A NUL-terminated, space-trimmed string of NBytes bytes at byte offset At.
    std::string Text(std::size_t At, std::size_t NBytes) const {
        if (At >= mSize)
            return {};
        const std::size_t n = std::min(NBytes, mSize - At);
        std::string s(mData + At, n);
        const auto nul = s.find('\0');
        if (nul != std::string::npos)
            s.resize(nul);
        const auto first = s.find_first_not_of(" \t\r\n\v\f");
        if (first == std::string::npos)
            return {};
        const auto last = s.find_last_not_of(" \t\r\n\v\f");
        return s.substr(first, last - first + 1);
    }
};

// (word size, swap) of a control block, trying 4/8-byte words in both orders.
std::optional<std::pair<int, bool>> d3_sniff(const char* pHead, std::size_t Size) {
    const bool little = std::endian::native == std::endian::little;
    for (const auto& [ws, file_little] :
         std::array<std::pair<int, bool>, 4>{{{4, true}, {4, false}, {8, true}, {8, false}}}) {
        if (Size < static_cast<std::size_t>(64 * ws))
            continue;
        const D3Words w{pHead, Size, ws, file_little != little};
        std::int64_t filetype = w.IntAt(11);
        if (filetype > 1000)
            filetype -= 1000;
        const std::int64_t ndim = w.IntAt(15);
        if ((filetype == 1 || filetype == 4 || filetype == 5 || filetype == 11) && ndim >= 2 &&
            ndim <= 9 && w.IntAt(16) >= 0)
            return std::make_pair(ws, file_little != little);
    }
    return std::nullopt;
}

// --- the family ----------------------------------------------------------------------

bool d3_all_digits(const std::string& rS) {
    return !rS.empty() &&
           std::all_of(rS.begin(), rS.end(), [](char c) { return c >= '0' && c <= '9'; });
}

// The base file followed by `<name>NN` members, by number (then name).
std::vector<std::string> d3_family(const std::string& rPath) {
    const fs::path path(rPath);
    const std::string base = path.filename().string();
    fs::path folder = path.parent_path();
    std::vector<std::pair<std::pair<unsigned long long, std::string>, std::string>> found;
    std::error_code ec;
    for (fs::directory_iterator it(folder.empty() ? fs::path(".") : folder, ec), end;
         !ec && it != end; it.increment(ec)) {
        const std::string name = it->path().filename().string();
        if (name.size() <= base.size() || name.compare(0, base.size(), base) != 0)
            continue;
        const std::string suffix = name.substr(base.size());
        if (!d3_all_digits(suffix))
            continue;
        std::error_code ec2;
        if (!fs::is_regular_file(it->path(), ec2))
            continue;
        unsigned long long number = 0;
        for (char c : suffix)
            number = number * 10 + static_cast<unsigned long long>(c - '0');
        found.push_back({{number, name}, (folder / name).string()});
    }
    std::sort(found.begin(), found.end());
    std::vector<std::string> out{rPath};
    for (auto& f : found)
        out.push_back(std::move(f.second));
    return out;
}

std::string d3_read_bytes(const std::string& rPath, std::size_t Offset, std::size_t N) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        d3_fail("cannot open '" + rPath + "'");
    in.seekg(static_cast<std::streamoff>(Offset));
    std::string out(N, '\0');
    in.read(out.data(), static_cast<std::streamsize>(N));
    out.resize(static_cast<std::size_t>(in.gcount()));
    return out;
}

// One past the last non-zero byte of a file (0 when all bytes are zero).
std::size_t d3_last_nonzero(const std::string& rPath, std::size_t Size) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        d3_fail("cannot open '" + rPath + "'");
    constexpr std::size_t chunk = std::size_t{1} << 16;
    std::string block;
    std::size_t end = Size;
    while (end > 0) {
        const std::size_t start = end > chunk ? end - chunk : 0;
        block.assign(end - start, '\0');
        in.seekg(static_cast<std::streamoff>(start));
        in.read(block.data(), static_cast<std::streamsize>(block.size()));
        for (std::size_t i = block.size(); i-- > 0;)
            if (block[i] != '\0')
                return start + i + 1;
        end = start;
    }
    return 0;
}

bool d3_file_sniffs(const std::string& rPath) {
    std::error_code ec;
    if (!fs::is_regular_file(rPath, ec))
        return false;
    const std::string head = d3_read_bytes(rPath, 0, 64 * 8);
    return d3_sniff(head.data(), head.size()).has_value();
}

// The base file when rPath is a numbered member of a family, else "".
std::string d3_continuation_base(const std::string& rPath) {
    const fs::path path(rPath);
    const std::string name = path.filename().string();
    std::size_t cut = name.size();
    while (cut > 0 && name[cut - 1] >= '0' && name[cut - 1] <= '9')
        --cut;
    if (cut == name.size() || cut == 0)
        return {};
    const fs::path base = path.parent_path() / name.substr(0, cut);
    std::error_code ec;
    return fs::is_regular_file(base, ec) ? base.string() : std::string();
}

// --- the control block ---------------------------------------------------------------

struct D3Header {
    std::map<std::string, std::int64_t> mRaw;
    std::size_t mBytes = 0;  // control block size
    int mFiletype = 1;
    bool mMaterialType = false, mRigidRoad = false, mRigidBodies = false,
         mReducedRigidBodies = false;
    std::int64_t mNodes = 0, mGlobals = 0;
    bool mMassScaling = false, mTemperature = false, mHeatFlux = false, mTemperatureLayers = false,
         mDisplacement = false, mVelocity = false, mAcceleration = false;
    std::int64_t mSolids = 0, mBeams = 0, mShells = 0, mTshells = 0;
    bool mSolidExtraNodes = false;
    std::int64_t mNv3d = 0, mNv1d = 0, mNv2d = 0, mNv3dt = 0, mNeiph = 0, mNeips = 0, mNeipb = 0,
                 mNt3d = 0;
    bool mElementDeletion = false, mNodeDeletion = false;
    std::int64_t mLayers = 0;
    bool mShellStress = false, mSolidStress = false, mShellPstrain = false, mSolidPstrain = false,
         mShellForces = false, mShellExtra = false;
    std::int64_t mParts = 0, mAirbags = 0, mAirbagSubver = 0, mShells8 = 0, mSph = 0;
    std::int64_t mSolids20 = 0, mSolids27 = 0;
    bool mQuadraticFull = false;  // QUADR > 0: a 27-node row lists all 27 nodes
    bool mTemperatureGradient = false, mResidualForces = false, mPlasticStrainTensor = false,
         mThermalStrainTensor = false, mElementStrain = false;

    std::int64_t R(const char* pName) const { return mRaw.at(pName); }

    std::int64_t SolidLayers() const {
        const std::int64_t base = 6 * mSolidStress + mSolidPstrain + mNeiph;
        return mNv3d / std::max<std::int64_t>(base, 1) >= 8 ? 8 : 1;
    }
};

D3Header d3_header(const D3Words& rW) {
    static const std::pair<int, const char*> words[] = {
        {10, "runtime"}, {11, "filetype"}, {12, "source_version"},
        {15, "ndim"},    {16, "numnp"},    {17, "icode"},
        {18, "nglbv"},   {19, "it"},       {20, "iu"},
        {21, "iv"},      {22, "ia"},       {23, "nel8"},
        {24, "nummat8"}, {25, "numds"},    {26, "numst"},
        {27, "nv3d"},    {28, "nel2"},     {29, "nummat2"},
        {30, "nv1d"},    {31, "nel4"},     {32, "nummat4"},
        {33, "nv2d"},    {34, "neiph"},    {35, "neips"},
        {36, "maxint"},  {37, "nmsph"},    {38, "ngpsph"},
        {39, "narbs"},   {40, "nelt"},     {41, "nummatt"},
        {42, "nv3dt"},   {43, "ioshl1"},   {44, "ioshl2"},
        {45, "ioshl3"},  {46, "ioshl4"},   {47, "ialemat"},
        {48, "ncfdv1"},  {49, "ncfdv2"},   {50, "nadapt"},
        {51, "nmmat"},   {52, "numfluid"}, {53, "inn"},
        {54, "npefg"},   {55, "nel48"},    {56, "idtdt"},
        {57, "extra"}};
    static const std::pair<int, const char*> extra_words[] = {
        {64, "nel20"},  {65, "nt3d"},    {66, "nel27"},  {67, "neipb"},
        {68, "nel21p"}, {69, "nel15t"},  {70, "soleng"}, {71, "nel20t"},
        {72, "nel40p"}, {73, "nel64"},   {74, "quadr"},  {75, "cubic"},
        {76, "tsheng"}, {77, "nbranch"}, {78, "penout"}, {79, "engout"}};
    D3Header h;
    for (const auto& [i, name] : words)
        h.mRaw[name] = rW.Int(static_cast<std::size_t>(i));
    const std::int64_t extra = h.mRaw["extra"];
    for (const auto& [i, name] : extra_words)
        h.mRaw[name] = (extra > 0 && i < 64 + extra) ? rW.Int(static_cast<std::size_t>(i)) : 0;
    h.mBytes = static_cast<std::size_t>(64 + std::max<std::int64_t>(extra, 0)) *
               static_cast<std::size_t>(rW.mWs);

    std::int64_t filetype = h.R("filetype");
    if (filetype > 1000)
        filetype -= 1000;
    if (filetype == 4)
        d3_fail("a intfor file is not read (only d3plot, d3part, d3eigv)");
    if (filetype != 1 && filetype != 5 && filetype != 11)
        d3_fail("unknown file type " + std::to_string(h.R("filetype")));
    h.mFiletype = static_cast<int>(filetype);
    if (h.R("nmmat") == kD3FemzipNmmat)
        d3_fail("the file is femzip-compressed; decompress it first");

    const std::int64_t ndim = h.R("ndim");
    h.mMaterialType = ndim == 5 || ndim == 7;
    h.mRigidRoad = ndim == 6 || ndim == 7 || ndim == 9;
    h.mRigidBodies = ndim == 8 || ndim == 9;
    h.mReducedRigidBodies = ndim == 9;
    if (ndim == 2)
        d3_fail("two-dimensional databases (NDIM = 2) are not read");
    h.mNodes = h.R("numnp");
    h.mGlobals = h.R("nglbv");

    const std::int64_t it = h.R("it");
    h.mMassScaling = d3_digit(it, 1) == 1;
    const std::int64_t it0 = d3_digit(it, 0);
    h.mTemperature = it0 >= 1 && it0 <= 3;
    h.mHeatFlux = it0 == 2 || it0 == 3;
    h.mTemperatureLayers = it0 == 3;
    h.mDisplacement = h.R("iu") != 0;
    h.mVelocity = h.R("iv") != 0;
    h.mAcceleration = h.R("ia") != 0;

    // NEL8 < 0 flags ten-node tetrahedra; its magnitude is the count (and
    // INT64_MIN, whose magnitude does not fit, is corrupt).
    const std::int64_t nel8 = h.R("nel8");
    if (nel8 == std::numeric_limits<std::int64_t>::min())
        d3_fail("the control block's NEL8 is corrupt");
    h.mSolids = nel8 < 0 ? -nel8 : nel8;
    h.mSolidExtraNodes = h.R("nel8") < 0;
    h.mBeams = h.R("nel2");
    h.mShells = h.R("nel4");
    h.mTshells = h.R("nelt");
    h.mNv3d = h.R("nv3d");
    h.mNv1d = h.R("nv1d");
    h.mNv2d = h.R("nv2d");
    h.mNv3dt = h.R("nv3dt");
    h.mNeiph = h.R("neiph");
    h.mNeips = h.R("neips");
    h.mNeipb = h.R("neipb");
    h.mNt3d = h.R("nt3d");

    const std::int64_t maxint = h.R("maxint");
    // Integration points per shell layer, encoded in MAXINT's magnitude; a
    // value whose magnitude does not fit is corrupt, not a layer count.
    if (maxint < -(std::int64_t{1} << 40) || maxint > (std::int64_t{1} << 40))
        d3_fail("the control block's MAXINT is corrupt");
    h.mElementDeletion = maxint <= -10000;
    h.mNodeDeletion = maxint > -10000 && maxint < 0;
    h.mLayers = maxint <= -10000 ? std::abs(maxint) - 10000 : std::abs(maxint);

    h.mShellStress = h.R("ioshl1") == 1000;
    h.mSolidStress = h.R("ioshl1") == 999 || h.R("ioshl1") == 1000;
    h.mShellPstrain = h.R("ioshl2") == 1000;
    h.mSolidPstrain = h.R("ioshl2") == 999 || h.R("ioshl2") == 1000;
    h.mShellForces = h.R("ioshl3") == 1000;
    h.mShellExtra = h.R("ioshl4") == 1000;

    if (h.R("ncfdv1") != 0)
        d3_fail("CFD or multi-solver data (NCFDV1 != 0) is not read");
    if (h.R("nadapt") != 0)
        d3_fail("adaptive-remeshing databases (NADAPT != 0) are not read");
    h.mParts = h.R("nmmat");
    const std::int64_t npefg = h.R("npefg");
    h.mAirbags = (npefg > 0 && npefg <= 10000000) ? npefg % 1000 : 0;
    h.mAirbagSubver = h.mAirbags ? npefg / 1000 : 0;
    h.mShells8 = h.R("nel48");
    h.mSph = h.R("nmsph");

    const std::int64_t idtdt = h.R("idtdt");
    h.mTemperatureGradient = d3_digit(idtdt, 0) == 1;
    h.mResidualForces = d3_digit(idtdt, 1) == 1;
    h.mPlasticStrainTensor = d3_digit(idtdt, 2) == 1;
    h.mThermalStrainTensor = d3_digit(idtdt, 3) == 1;
    // In doubles: the words are unchecked here, and int64 products of corrupt
    // ones overflow (d3_state_words rejects such a header afterwards).
    const double layer_vars = 6.0 * static_cast<double>(h.mShellStress) +
                              static_cast<double>(h.mShellPstrain) + static_cast<double>(h.mNeips);
    const double layers = static_cast<double>(h.mLayers);
    if (idtdt > 100)
        h.mElementStrain = d3_digit(idtdt, 4) == 1;
    else if (h.mNv2d > 0)
        h.mElementStrain = static_cast<double>(h.mNv2d) - layers * layer_vars -
                               8.0 * static_cast<double>(h.mShellForces) -
                               4.0 * static_cast<double>(h.mShellExtra) >
                           1.0;
    else if (h.mNv3dt > 0)
        h.mElementStrain = static_cast<double>(h.mNv3dt) - layers * layer_vars > 1.0;

    // 20- and 27-node hexahedra are read; the node order of the others
    // (21-node wedges, 15-node tetrahedra, the cubic solids) is not documented
    // in the manuals available.
    h.mSolids20 = h.R("nel20");
    h.mSolids27 = h.R("nel27");
    h.mQuadraticFull = h.R("quadr") > 0;
    for (const char* key : {"nel21p", "nel15t", "nel20t", "nel40p", "nel64"})
        if (h.R(key) > 0) {
            std::string upper = key;
            for (char& c : upper)
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            d3_fail("higher-order solids (" + upper + " = " + std::to_string(h.R(key)) +
                    ") are not read");
        }
    // Every node, element and part is at least a word of the base file; a
    // larger count is a corrupt control block, and each sizes an id table.
    for (const std::int64_t n : {h.mNodes, h.mSolids, h.mBeams, h.mShells, h.mTshells, h.mParts})
        if (n < 0 || static_cast<std::uint64_t>(n) > rW.NumWords())
            d3_fail("the control block counts more entities than the file holds");
    return h;
}

// --- the geometry --------------------------------------------------------------------

struct D3Geometry {
    std::int64_t mRigidShells = 0;
    std::vector<std::int64_t> mPartMattype;
    std::int64_t mSphVars = 0;
    bool mHasAirbag = false;
    std::int64_t mAirbagGeom = 0, mAirbagVar = 0, mAirbagParticles = 0, mAirbagStateGeom = 0;
    std::vector<double> mCoords;
    std::vector<std::int64_t> mSolids, mTshells, mBeams, mShells;  // 9/9/6/5 words per element
    std::vector<std::int64_t> mSolidExtra;                         // 2 per solid
    std::vector<std::int64_t> mNodeIds, mSolidIds, mBeamIds, mShellIds, mTshellIds;
    std::vector<std::int64_t> mPartIds;
    bool mHasPartIds = false;
    std::int64_t mRigidBodies = 0;       // from the numbering section
    std::int64_t mRigidBodyMotions = 0;  // rigid bodies with motion in the states
    std::int64_t mRoads = 0;
    std::vector<std::int64_t> mShell8;           // 5 words per 8-node shell
    std::vector<std::int64_t> mSphFlags;         // ISPHFG(1..10) and the history count
    std::vector<std::int64_t> mSph;              // (node, material) per particle, 1-based
    std::vector<std::int64_t> mAirbagTypes;      // 1 integer, 2 real, per variable
    std::vector<std::string> mAirbagNames;       // geometry, particle, then bag variables
    std::vector<std::int64_t> mAirbagGeomData;   // per bag: first particle, count, id...
    std::vector<std::int64_t> mRigidBodyParts;   // per rigid body its part, 1-based
    std::vector<std::int64_t> mRoadNodeIds;      // rigid road nodes' ids
    std::vector<double> mRoadCoords;             // and coordinates
    std::vector<std::int64_t> mRoadSegments;     // 4 road node ids per segment
    std::vector<std::int64_t> mRoadSegmentRoad;  // each segment's road id
    std::vector<std::int64_t> mSolid20;          // 13 words per 20-node hexahedron
    std::vector<std::int64_t> mSolid27;          // 20 (or 28) words per 27-node one
    std::size_t mSolid27Width = 20;
    std::map<std::int64_t, std::string> mPartTitles;
    std::vector<std::int64_t> mPartTitleOrder;
    std::size_t mEnd = 0;  // words
};

std::vector<std::int64_t> d3_iota(std::int64_t N) {
    std::vector<std::int64_t> v(static_cast<std::size_t>(std::max<std::int64_t>(N, 0)));
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = static_cast<std::int64_t>(i) + 1;
    return v;
}

D3Geometry d3_geometry(const D3Words& rW, const D3Header& rH) {
    D3Geometry g;
    const auto sz = [](std::int64_t N) {
        return static_cast<std::size_t>(std::max<std::int64_t>(N, 0));
    };
    std::size_t pos = rH.mBytes / static_cast<std::size_t>(rW.mWs);
    const auto ws = static_cast<std::size_t>(rW.mWs);
    if (rH.mMaterialType) {
        g.mRigidShells = rW.Int(pos);
        const std::int64_t nummat = rW.Int(pos + 1);
        if (nummat != rH.mParts)
            d3_fail("the material type section lists " + std::to_string(nummat) +
                    " parts, the header " + std::to_string(rH.mParts));
        g.mPartMattype = rW.Ints(pos + 2, sz(rH.mParts));
        pos += 2 + sz(rH.mParts);
    }
    if (rH.R("ialemat") > 0)
        pos += sz(rH.R("ialemat"));
    if (rH.mSph > 0) {
        const auto flags = rW.Ints(pos, 11);
        for (const std::int64_t f : flags)  // counts of SPH variables: small
            if (f < -(1 << 20) || f > (1 << 20))
                d3_fail("the SPH flags hold an implausible count");
        // ISPHFG(1) = 10 (newer releases) leaves ISPHFG(11) undefined
        const std::int64_t history = flags[0] == 10 ? 0 : flags[10];
        std::int64_t sum = 0;
        for (int i = 1; i < 8; ++i)
            sum += flags[static_cast<std::size_t>(i)];
        g.mSphVars = sum + std::abs(flags[8]) + flags[9] + history + 1;
        g.mSphFlags.assign(flags.begin(), flags.begin() + 10);
        g.mSphFlags.push_back(history);
        pos += sz(flags[0]);
    }
    if (rH.mAirbags) {
        const auto head = rW.Ints(pos, 4);
        g.mAirbagGeom = head[0];
        g.mAirbagVar = head[1];
        g.mAirbagParticles = head[2];
        g.mAirbagStateGeom = head[3];
        // Each sizes an allocation below (particle points, geometry-data
        // reads); a value the file cannot back is corrupt, matching the
        // control block's own node/element/part counts (above).
        for (const std::int64_t v :
             {g.mAirbagGeom, g.mAirbagVar, g.mAirbagParticles, g.mAirbagStateGeom})
            if (v < 0 || static_cast<std::uint64_t>(v) > rW.NumWords())
                d3_fail("the airbag header counts more entities than the file holds");
        g.mHasAirbag = true;
        pos += 4;
        if (rH.mAirbagSubver == 4)
            pos += 1;
        // type codes (1 integer, 2 real), then 8-word names, one character
        // per word: geometry, particle, then bag variables
        const std::size_t nvars = sz(g.mAirbagGeom + g.mAirbagVar + g.mAirbagStateGeom);
        g.mAirbagTypes = rW.Ints(pos, nvars);
        pos += nvars;
        for (std::size_t v = 0; v < nvars; ++v) {
            std::string name;
            for (std::int64_t c : rW.Ints(pos + 8 * v, 8))
                name.push_back(static_cast<char>(c & 0xFF));
            const std::string blank(" \0", 2);
            const std::size_t a = name.find_first_not_of(blank);
            const std::size_t b = name.find_last_not_of(blank);
            g.mAirbagNames.push_back(a == std::string::npos ? "" : name.substr(a, b - a + 1));
        }
        pos += 8 * nvars;
    }

    const std::size_t n = sz(rH.mNodes);
    g.mCoords = rW.Floats(pos, 3 * n);
    pos += 3 * n;
    g.mSolids = rW.Ints(pos, 9 * sz(rH.mSolids));
    pos += 9 * sz(rH.mSolids);
    if (rH.mSolidExtraNodes) {
        g.mSolidExtra = rW.Ints(pos, 2 * sz(rH.mSolids));
        pos += 2 * sz(rH.mSolids);
    }
    g.mTshells = rW.Ints(pos, 9 * sz(rH.mTshells));
    pos += 9 * sz(rH.mTshells);
    g.mBeams = rW.Ints(pos, 6 * sz(rH.mBeams));
    pos += 6 * sz(rH.mBeams);
    g.mShells = rW.Ints(pos, 5 * sz(rH.mShells));
    pos += 5 * sz(rH.mShells);

    g.mNodeIds = d3_iota(rH.mNodes);
    g.mSolidIds = d3_iota(rH.mSolids);
    g.mBeamIds = d3_iota(rH.mBeams);
    g.mShellIds = d3_iota(rH.mShells);
    g.mTshellIds = d3_iota(rH.mTshells);
    const std::int64_t narbs = rH.R("narbs");
    if (narbs > 0) {
        const std::size_t start = pos;
        const auto head = rW.Ints(pos, 10);
        pos += 10;
        const std::int64_t nsort = head[0];
        std::int64_t nparts_ids = rH.mParts;
        if (nsort < 0) {
            const auto extra = rW.Ints(pos, 6);
            pos += 6;
            nparts_ids = extra[3];
            g.mRigidBodies = extra[4];
        }
        std::vector<std::vector<std::int64_t>> ids;
        for (int k = 5; k < 10; ++k) {
            const std::size_t count = sz(head[static_cast<std::size_t>(k)]);
            ids.push_back(rW.Ints(pos, count));
            pos += count;
        }
        if (ids[0].size() == n)
            g.mNodeIds = ids[0];
        if (ids[1].size() == sz(rH.mSolids))
            g.mSolidIds = ids[1];
        if (ids[2].size() == sz(rH.mBeams))
            g.mBeamIds = ids[2];
        if (ids[3].size() == sz(rH.mShells))
            g.mShellIds = ids[3];
        if (ids[4].size() == sz(rH.mTshells))
            g.mTshellIds = ids[4];
        if (nsort < 0 && nparts_ids == rH.mParts) {
            g.mPartIds = rW.Ints(pos, sz(rH.mParts));
            g.mHasPartIds = true;
        }
        pos = start + sz(narbs);
    }

    if (rH.mRigidBodies) {
        const std::int64_t nrigid = rW.Int(pos);
        pos += 1;
        for (std::int64_t r = 0; r < nrigid; ++r) {
            g.mRigidBodyParts.push_back(rW.Int(pos));
            const std::int64_t numnodr = rW.Int(pos + 1);
            pos += 2 + sz(numnodr);
            const std::int64_t numnoda = rW.Int(pos);
            pos += 1 + sz(numnoda);
        }
        g.mRigidBodyMotions = nrigid;
    }
    // SPH particles: (node, material), 1-based.
    if (rH.mSph > 0) {
        g.mSph = rW.Ints(pos, 2 * sz(rH.mSph));
        pos += 2 * sz(rH.mSph);
    }
    // Airbags: per bag first particle, particle count, id, gas mixtures (and
    // chambers).
    if (g.mHasAirbag) {
        g.mAirbagGeomData = rW.Ints(pos, sz(rH.mAirbags * g.mAirbagGeom));
        pos += sz(rH.mAirbags * g.mAirbagGeom);
    }
    // Rigid roads: their own nodes (ids, coordinates) and per surface its id
    // and its 4-node segments (road node ids).
    if (rH.mRigidRoad) {
        const auto head = rW.Ints(pos, 3);
        const std::size_t nnode = sz(head[0]);
        pos += 4;
        g.mRoadNodeIds = rW.Ints(pos, nnode);
        pos += nnode;
        g.mRoadCoords = rW.Floats(pos, 3 * nnode);
        pos += 3 * nnode;
        for (std::int64_t s = 0; s < head[2]; ++s) {
            const std::int64_t road = rW.Int(pos);
            const std::size_t nseg = sz(rW.Int(pos + 1));
            pos += 2;
            const auto segs = rW.Ints(pos, 4 * nseg);
            g.mRoadSegments.insert(g.mRoadSegments.end(), segs.begin(), segs.end());
            g.mRoadSegmentRoad.insert(g.mRoadSegmentRoad.end(), nseg, road);
            pos += 4 * nseg;
        }
        g.mRoads = head[2];
    }

    // The ten-node solids' extra nodes follow the 8-node connectivity (database
    // manual); lasso-python also writes and reads a second, identical copy
    // here, which is skipped when present.
    if (rH.mSolidExtraNodes && rH.mSolids > 0) {
        const std::size_t n2 = 2 * sz(rH.mSolids);
        if ((pos + n2) * ws <= rW.mSize && rW.Ints(pos, n2) == g.mSolidExtra)
            pos += n2;
    }
    if (rH.mShells8 > 0) {
        g.mShell8 = rW.Ints(pos, 5 * sz(rH.mShells8));
        pos += 5 * sz(rH.mShells8);
    }
    // 20-node hexahedra: the solid's index and its 12 edge nodes; 27-node
    // ones: the index and nodes 9-27 (all 27 when QUADR > 0), in LS-DYNA's
    // order, which is VTK's (keyword manual, *ELEMENT_SOLID).
    if (rH.mSolids20 > 0) {
        g.mSolid20 = rW.Ints(pos, 13 * sz(rH.mSolids20));
        pos += 13 * sz(rH.mSolids20);
    }
    g.mSolid27Width = rH.mQuadraticFull ? 28 : 20;
    if (rH.mSolids27 > 0) {
        g.mSolid27 = rW.Ints(pos, g.mSolid27Width * sz(rH.mSolids27));
        pos += g.mSolid27Width * sz(rH.mSolids27);
    }

    if (pos < rW.NumWords() && rW.FloatAt(pos) == kD3EofMarker) {
        pos += 1;
        while (pos < rW.NumWords()) {
            const std::int64_t ntype = rW.Int(pos);
            if (ntype == 90000) {
                pos += 1 + 72 / ws;
            } else if (ntype == 90001 || ntype == 90002 || ntype == 90020) {
                const std::int64_t count = rW.Int(pos + 1);
                pos += 2;
                const std::size_t entry = ws + 72;
                for (std::int64_t k = 0; k < count; ++k) {
                    const std::size_t at = pos * ws + static_cast<std::size_t>(k) * entry;
                    if (ntype == 90001) {
                        const std::int64_t pid = rW.Int(at / ws);
                        if (g.mPartTitles.emplace(pid, rW.Text(at + ws, 72)).second)
                            g.mPartTitleOrder.push_back(pid);
                    }
                }
                pos += (sz(count) * entry) / ws;
            } else {
                break;
            }
            if (pos < rW.NumWords() && rW.FloatAt(pos) == kD3EofMarker)
                pos += 1;
        }
    }
    g.mEnd = pos;
    return g;
}

// --- the state layout ----------------------------------------------------------------

std::vector<std::pair<std::string, std::size_t>> d3_node_vars(const D3Header& rH) {
    std::vector<std::pair<std::string, std::size_t>> out;
    if (rH.mDisplacement)
        out.emplace_back("coordinates", 3);
    if (rH.mTemperature)
        out.emplace_back("temperature", rH.mTemperatureLayers ? 3 : 1);
    if (rH.mHeatFlux)
        out.emplace_back("heat_flux", 3);
    if (rH.mMassScaling)
        out.emplace_back("mass_scaling", 1);
    if (rH.mTemperatureGradient)
        out.emplace_back("temperature_gradient", 1);
    if (rH.mResidualForces) {
        out.emplace_back("residual_forces", 3);
        out.emplace_back("residual_moments", 3);
    }
    if (rH.mVelocity)
        out.emplace_back("velocity", 3);
    if (rH.mAcceleration)
        out.emplace_back("acceleration", 3);
    return out;
}

std::int64_t d3_state_words(const D3Header& rH, const D3Geometry& rG) {
    // Every term is a count from the control block. A negative one, or a sum
    // past int64, would let a state's element blocks index outside the state
    // record the size is checked against -- so both are corruption.
    constexpr std::int64_t kMax = std::numeric_limits<std::int64_t>::max();
    const auto nonneg = [](std::int64_t V) {
        if (V < 0)
            d3_fail("the control block holds a negative count");
        return V;
    };
    const auto mul = [&](std::int64_t A, std::int64_t B) {
        nonneg(A);
        nonneg(B);
        if (A != 0 && B > kMax / A)
            d3_fail("the state size overflows");
        return A * B;
    };
    std::int64_t n = 0;
    const auto add = [&](std::int64_t V) {
        if (nonneg(V) > kMax - n)
            d3_fail("the state size overflows");
        n += V;
    };
    add(1);
    add(rH.mGlobals);
    std::int64_t comps = 0;
    for (const auto& v : d3_node_vars(rH))
        comps += static_cast<std::int64_t>(v.second);
    add(mul(comps, rH.mNodes));
    add(mul(rH.mNt3d, rH.mSolids));
    add(mul(rH.mSolids, rH.mNv3d));
    add(mul(rH.mTshells, rH.mNv3dt));
    add(mul(rH.mBeams, rH.mNv1d));
    if (rG.mRigidShells > rH.mShells)
        d3_fail("more rigid shells than shells");
    add(mul(rH.mShells - rG.mRigidShells, rH.mNv2d));
    add(mul(rH.mSph, rG.mSphVars));
    if (rH.mNodeDeletion)
        add(rH.mNodes);
    else if (rH.mElementDeletion) {
        add(rH.mBeams);
        add(rH.mShells);
        add(rH.mSolids);
        add(rH.mTshells);
    }
    if (rG.mHasAirbag) {
        add(mul(rH.mAirbags, rG.mAirbagStateGeom));
        add(mul(rG.mAirbagParticles, rG.mAirbagVar));
    }
    add(mul(rG.mRoads, 6));
    if (rH.mRigidBodies)
        add(mul(rG.mRigidBodyMotions, rH.mReducedRigidBodies ? 12 : 24));
    return n;
}

// The family: header, geometry, and where each state lies.
struct D3File {
    std::string mBase;  // the base file's bytes
    D3Words mWords;
    D3Header mHeader;
    D3Geometry mGeometry;
    std::size_t mStateWords = 0;
    std::vector<std::pair<std::string, std::size_t>> mStates;  // (file, byte offset)

    explicit D3File(const std::string& rPath) {
        std::error_code ec;
        if (!fs::is_regular_file(rPath, ec))
            d3_fail("'" + rPath + "' does not exist");
        const std::string base = d3_continuation_base(rPath);
        if (!base.empty() && d3_file_sniffs(base))
            d3_fail("'" + fs::path(rPath).filename().string() + "' continues the family of '" +
                    fs::path(base).filename().string() + "'; open the base file");
        const std::size_t size = static_cast<std::size_t>(fs::file_size(rPath, ec));
        mBase = d3_read_bytes(rPath, 0, size);
        std::string lower = rPath;
        for (char& c : lower)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (lower.size() >= 3 && lower.compare(lower.size() - 3, 3, ".fz") == 0)
            d3_fail("the file is femzip-compressed; decompress it first");
        const auto sniffed = d3_sniff(mBase.data(), mBase.size());
        if (!sniffed)
            d3_fail("'" + rPath + "' is not a d3plot file (no plausible control block)");
        mWords = D3Words{mBase.data(), mBase.size(), sniffed->first, sniffed->second};
        mHeader = d3_header(mWords);
        mGeometry = d3_geometry(mWords, mHeader);
        const std::int64_t words = d3_state_words(mHeader, mGeometry);
        if (words <= 0)
            d3_fail("the header describes an empty state");
        mStateWords = static_cast<std::size_t>(words);
        const std::size_t state_bytes = mStateWords * static_cast<std::size_t>(mWords.mWs);

        const auto members = d3_family(rPath);
        for (std::size_t k = 0; k < members.size(); ++k) {
            const std::size_t fsize = static_cast<std::size_t>(fs::file_size(members[k], ec));
            const std::size_t start =
                k == 0 ? mGeometry.mEnd * static_cast<std::size_t>(mWords.mWs) : 0;
            // the last non-zero *word*: a big-endian value can end in zero bytes
            const std::size_t ws = static_cast<std::size_t>(mWords.mWs);
            const std::size_t last = (d3_last_nonzero(members[k], fsize) + ws - 1) / ws * ws;
            std::size_t count = last > start ? (last - start) / state_bytes : 0;
            if (mHeader.mFiletype == 11 && k > 0 && count == 0 && state_bytes <= fsize)
                count = 1;
            for (std::size_t i = 0; i < count; ++i)
                mStates.emplace_back(members[k], start + i * state_bytes);
        }
    }

    std::vector<double> State(std::size_t Index) const {
        const std::size_t bytes = mStateWords * static_cast<std::size_t>(mWords.mWs);
        const std::string raw = d3_read_bytes(mStates[Index].first, mStates[Index].second, bytes);
        if (raw.size() < bytes)
            d3_fail("state " + std::to_string(Index) + " is truncated");
        const D3Words w{raw.data(), raw.size(), mWords.mWs, mWords.mSwap};
        return w.Floats(0, mStateWords);
    }

    // The words of state Index read as integers.
    std::vector<std::int64_t> StateInts(std::size_t Index) const {
        const std::size_t bytes = mStateWords * static_cast<std::size_t>(mWords.mWs);
        const std::string raw = d3_read_bytes(mStates[Index].first, mStates[Index].second, bytes);
        if (raw.size() < bytes)
            d3_fail("state " + std::to_string(Index) + " is truncated");
        const D3Words w{raw.data(), raw.size(), mWords.mWs, mWords.mSwap};
        return w.Ints(0, mStateWords);
    }

    std::vector<double> Times() const {
        std::vector<double> out;
        out.reserve(mStates.size());
        for (const auto& [member, offset] : mStates) {
            const std::string raw =
                d3_read_bytes(member, offset, static_cast<std::size_t>(mWords.mWs));
            if (raw.size() < static_cast<std::size_t>(mWords.mWs))
                d3_fail("a state is truncated");
            const D3Words w{raw.data(), raw.size(), mWords.mWs, mWords.mSwap};
            out.push_back(w.FloatAt(0));
        }
        return out;
    }
};

// --- the mesh ------------------------------------------------------------------------

enum D3Family : int {
    kD3Solid = 0,
    kD3Tshell = 1,
    kD3Beam = 2,
    kD3Shell = 3,
    kD3Sph = 4,
    kD3Airbag = 5,
    kD3Road = 6
};
constexpr std::size_t kD3Families = 7;

struct D3Block {
    std::string mType;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;
    int mFamily = 0;
    std::vector<std::size_t> mElems;
    std::vector<std::int64_t> mParts;  // 0-based part indices
};

struct D3Cells {
    std::vector<D3Block> mBlocks;
    // per family: (block, row) of each element
    std::array<std::vector<std::pair<std::size_t, std::size_t>>, kD3Families> mWhere;
};

int d3_dim(const std::string& rType) {
    if (rType == "vertex")
        return 0;
    if (rType == "line")
        return 1;
    if (rType == "triangle" || rType == "quad" || rType == "quad8")
        return 2;
    return 3;
}

D3Cells d3_cells(const D3Header& rH, const D3Geometry& rG) {
    D3Cells out;
    const auto add = [&](int Family, const std::vector<std::int64_t>& rRows, std::size_t Width,
                         std::size_t PartWord, const auto& rMake) {
        const std::size_t count = Width ? rRows.size() / Width : 0;
        std::map<std::string, std::size_t> index;
        auto& where = out.mWhere[static_cast<std::size_t>(Family)];
        where.resize(count);
        for (std::size_t e = 0; e < count; ++e) {
            const std::int64_t* row = rRows.data() + e * Width;
            auto [type, nodes] = rMake(e, row);
            auto it = index.find(type);
            if (it == index.end()) {
                it = index.emplace(type, out.mBlocks.size()).first;
                D3Block b;
                b.mType = type;
                b.mNodes = nodes.size();
                b.mFamily = Family;
                out.mBlocks.push_back(std::move(b));
            }
            D3Block& b = out.mBlocks[it->second];
            where[e] = {it->second, b.mElems.size()};
            b.mConn.insert(b.mConn.end(), nodes.begin(), nodes.end());
            b.mElems.push_back(e);
            b.mParts.push_back(d3_zero_based(row[PartWord]));
        }
    };
    // 20- and 27-node hexahedra: the corners, then nodes 9.. in VTK's order.
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> solid20, solid27;
    for (std::size_t i = 0; i + 13 <= rG.mSolid20.size(); i += 13) {
        auto& v = solid20[d3_zero_based(rG.mSolid20[i])];
        for (std::size_t k = 1; k < 13; ++k)
            v.push_back(d3_zero_based(rG.mSolid20[i + k]));
    }
    const std::size_t w27 = rG.mSolid27Width;
    for (std::size_t i = 0; i + w27 <= rG.mSolid27.size(); i += w27) {
        auto& v = solid27[d3_zero_based(rG.mSolid27[i])];
        for (std::size_t k = 1; k < w27; ++k)
            v.push_back(d3_zero_based(rG.mSolid27[i + k]));
    }
    const auto solid = [&](std::size_t E, const std::int64_t* pRow) {
        std::array<std::int64_t, 8> nodes;
        for (std::size_t k = 0; k < 8; ++k)
            nodes[k] = d3_zero_based(pRow[k]);
        const auto e = static_cast<std::int64_t>(E);
        if (const auto it = solid20.find(e); it != solid20.end()) {
            std::vector<std::int64_t> v(nodes.begin(), nodes.end());
            v.insert(v.end(), it->second.begin(), it->second.end());
            return std::make_pair(std::string("hexahedron20"), v);
        }
        if (const auto it = solid27.find(e); it != solid27.end()) {
            if (it->second.size() == 27)
                return std::make_pair(std::string("hexahedron27"), it->second);
            std::vector<std::int64_t> v(nodes.begin(), nodes.end());
            v.insert(v.end(), it->second.begin(), it->second.end());
            return std::make_pair(std::string("hexahedron27"), v);
        }
        if (!rG.mSolidExtra.empty() && rG.mSolidExtra[2 * E] > 0 && rG.mSolidExtra[2 * E + 1] > 0) {
            std::vector<std::int64_t> v(nodes.begin(), nodes.end());
            v.push_back(d3_zero_based(rG.mSolidExtra[2 * E]));
            v.push_back(d3_zero_based(rG.mSolidExtra[2 * E + 1]));
            return std::make_pair(std::string("tetra10"), v);
        }
        const detail::CollapsedBrick c = detail::collapse_brick(nodes);
        return std::make_pair(std::string(c.mType), c.mNodes);
    };
    const auto tshell = [&](std::size_t /*E*/, const std::int64_t* pRow) {
        std::array<std::int64_t, 8> nodes;
        for (std::size_t k = 0; k < 8; ++k)
            nodes[k] = d3_zero_based(pRow[k]);
        const detail::CollapsedBrick c = detail::collapse_brick(nodes);
        return std::make_pair(std::string(c.mType), c.mNodes);
    };
    const auto beam = [](std::size_t /*E*/, const std::int64_t* pRow) {
        return std::make_pair(
            std::string("line"),
            std::vector<std::int64_t>{d3_zero_based(pRow[0]), d3_zero_based(pRow[1])});
    };
    std::unordered_map<std::int64_t, std::array<std::int64_t, 4>> shell8;
    for (std::size_t i = 0; i + 5 <= rG.mShell8.size(); i += 5)
        shell8[d3_zero_based(rG.mShell8[i])] = {
            d3_zero_based(rG.mShell8[i + 1]), d3_zero_based(rG.mShell8[i + 2]),
            d3_zero_based(rG.mShell8[i + 3]), d3_zero_based(rG.mShell8[i + 4])};
    const auto shell = [&](std::size_t E, const std::int64_t* pRow) {
        std::vector<std::int64_t> nodes{d3_zero_based(pRow[0]), d3_zero_based(pRow[1]),
                                        d3_zero_based(pRow[2]), d3_zero_based(pRow[3])};
        const auto it = shell8.find(static_cast<std::int64_t>(E));
        if (it != shell8.end()) {
            nodes.insert(nodes.end(), it->second.begin(), it->second.end());
            return std::make_pair(std::string("quad8"), nodes);
        }
        if (nodes[3] == nodes[2] || nodes[3] < 0) {
            nodes.resize(3);
            return std::make_pair(std::string("triangle"), nodes);
        }
        return std::make_pair(std::string("quad"), nodes);
    };
    add(kD3Solid, rG.mSolids, 9, 8, solid);
    add(kD3Tshell, rG.mTshells, 9, 8, tshell);
    add(kD3Beam, rG.mBeams, 6, 5, beam);
    add(kD3Shell, rG.mShells, 5, 4, shell);
    // SPH particles: vertices on their nodes, the material as the part.
    add(kD3Sph, rG.mSph, 2, 1, [](std::size_t /*E*/, const std::int64_t* pRow) {
        return std::make_pair(std::string("vertex"),
                              std::vector<std::int64_t>{d3_zero_based(pRow[0])});
    });
    // Airbag particles and rigid road nodes are points after the nodes.
    const auto nparticles = static_cast<std::size_t>(
        std::max<std::int64_t>(rG.mHasAirbag ? rG.mAirbagParticles : 0, 0));
    if (nparticles) {
        std::vector<std::int64_t> owner(nparticles, 0);
        const auto ngeom = static_cast<std::size_t>(std::max<std::int64_t>(rG.mAirbagGeom, 1));
        for (std::size_t b = 0; b * ngeom + 1 < rG.mAirbagGeomData.size(); ++b) {
            const std::int64_t first = rG.mAirbagGeomData[b * ngeom];
            const std::int64_t count = rG.mAirbagGeomData[b * ngeom + 1];
            if (count == 0)
                continue;
            if (first < 1 || static_cast<std::uint64_t>(first - 1) > nparticles || count < 0 ||
                static_cast<std::uint64_t>(count) >
                    nparticles - static_cast<std::size_t>(first - 1))
                d3_fail("an airbag references particles outside the particle table");
            const auto start = static_cast<std::size_t>(first - 1);
            for (std::size_t i = start; i < start + static_cast<std::size_t>(count); ++i)
                owner[i] = static_cast<std::int64_t>(b) + 1;
        }
        std::vector<std::int64_t> rows;
        for (std::size_t i = 0; i < nparticles; ++i) {
            rows.push_back(rH.mNodes + static_cast<std::int64_t>(i));
            rows.push_back(owner[i]);
        }
        add(kD3Airbag, rows, 2, 1, [](std::size_t /*E*/, const std::int64_t* pRow) {
            return std::make_pair(std::string("vertex"), std::vector<std::int64_t>{pRow[0]});
        });
    }
    if (!rG.mRoadSegments.empty()) {
        const std::int64_t base = rH.mNodes + static_cast<std::int64_t>(nparticles);
        std::unordered_map<std::int64_t, std::int64_t> index;
        for (std::size_t i = 0; i < rG.mRoadNodeIds.size(); ++i)
            index.emplace(rG.mRoadNodeIds[i], base + static_cast<std::int64_t>(i));
        std::vector<std::int64_t> rows;
        for (std::size_t s = 0; s < rG.mRoadSegmentRoad.size(); ++s) {
            for (std::size_t k = 0; k < 4; ++k) {
                const auto it = index.find(rG.mRoadSegments[4 * s + k]);
                if (it == index.end())
                    d3_fail("a rigid road segment names a node the road does not have");
                rows.push_back(it->second);
            }
            rows.push_back(rG.mRoadSegmentRoad[s]);
        }
        add(kD3Road, rows, 5, 4, [](std::size_t /*E*/, const std::int64_t* pRow) {
            return std::make_pair(std::string("quad"),
                                  std::vector<std::int64_t>{pRow[0], pRow[1], pRow[2], pRow[3]});
        });
    }
    return out;
}

std::vector<std::int64_t> d3_part_user_ids(const D3Header& rH, const D3Geometry& rG) {
    if (rG.mHasPartIds)
        return rG.mPartIds;
    if (static_cast<std::int64_t>(rG.mPartTitleOrder.size()) == rH.mParts)
        return rG.mPartTitleOrder;
    return d3_iota(rH.mParts);
}

NDArray d3_int_array(const std::vector<std::int64_t>& rValues) {
    NDArray a(DType::Int64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), a.As<std::int64_t>());
    return a;
}

NDArray d3_scalar(double V, DType T) {
    NDArray a(T, {1});
    if (T == DType::Float64)
        a.As<double>()[0] = V;
    else
        a.As<std::int64_t>()[0] = static_cast<std::int64_t>(V);
    return a;
}

Mesh d3_build_mesh(const D3File& rF, D3Cells& rCells) {
    const D3Header& h = rF.mHeader;
    const D3Geometry& g = rF.mGeometry;
    const auto nnodes = static_cast<std::size_t>(h.mNodes);
    // the nodes, then the airbag particles (placed by each state) and the
    // rigid road nodes
    const auto nparticles =
        static_cast<std::size_t>(std::max<std::int64_t>(g.mHasAirbag ? g.mAirbagParticles : 0, 0));
    const std::size_t nroad = g.mRoadNodeIds.size();
    const std::size_t npts = nnodes + nparticles + nroad;
    Mesh mesh;
    NDArray points(DType::Float64, {npts, 3});
    double* pp = points.As<double>();
    std::copy(g.mCoords.begin(), g.mCoords.end(), pp);
    std::fill(pp + 3 * nnodes, pp + 3 * (nnodes + nparticles), kD3Nan);
    std::copy(g.mRoadCoords.begin(), g.mRoadCoords.end(), pp + 3 * (nnodes + nparticles));
    mesh.AssignPoints(std::move(points));
    rCells = d3_cells(h, g);
    for (const D3Block& b : rCells.mBlocks)
        for (std::int64_t v : b.mConn)
            if (v < 0 || static_cast<std::size_t>(v) >= npts)
                d3_fail("an element references a node outside the node table");
    for (const D3Block& b : rCells.mBlocks) {
        NDArray conn(DType::Int64, {b.mElems.size(), b.mNodes});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
    }
    {
        std::vector<std::int64_t> nid(npts, -1);
        std::copy(g.mNodeIds.begin(), g.mNodeIds.end(), nid.begin());
        mesh.AddPointData("lsdyna:nid", d3_int_array(nid));
    }
    const auto part_ids = d3_part_user_ids(h, g);
    if (!g.mRigidBodyParts.empty()) {  // each rigid body's part (user id)
        std::vector<std::int64_t> rigid;
        rigid.reserve(g.mRigidBodyParts.size());
        for (std::int64_t p : g.mRigidBodyParts)
            rigid.push_back(p > 0 && static_cast<std::size_t>(p) <= part_ids.size()
                                ? part_ids[static_cast<std::size_t>(p - 1)]
                                : p);
        mesh.AddFieldData("lsdyna:rigid_body_part", d3_int_array(rigid));
    }
    if (rCells.mBlocks.empty())
        return mesh;

    // an SPH particle is named by its node, an airbag particle and a road
    // segment by their number
    std::vector<std::int64_t> sph_ids;
    for (std::size_t i = 0; i + 1 < g.mSph.size(); i += 2) {
        const std::int64_t node = d3_zero_based(g.mSph[i]);
        if (node < 0 || static_cast<std::size_t>(node) >= g.mNodeIds.size())
            d3_fail("an SPH particle references a node outside the node table");
        sph_ids.push_back(g.mNodeIds[static_cast<std::size_t>(node)]);
    }
    const std::vector<std::int64_t> particle_ids = d3_iota(static_cast<std::int64_t>(nparticles));
    const std::vector<std::int64_t> segment_ids =
        d3_iota(static_cast<std::int64_t>(g.mRoadSegmentRoad.size()));
    const std::array<const std::vector<std::int64_t>*, kD3Families> family_ids{
        &g.mSolidIds, &g.mTshellIds, &g.mBeamIds, &g.mShellIds,
        &sph_ids,     &particle_ids, &segment_ids};
    std::vector<std::int64_t> bag_ids;
    if (g.mAirbagGeom > 0)
        for (std::size_t b = 0;
             b * static_cast<std::size_t>(g.mAirbagGeom) + 2 < g.mAirbagGeomData.size(); ++b)
            bag_ids.push_back(g.mAirbagGeomData[b * static_cast<std::size_t>(g.mAirbagGeom) + 2]);
    // airbag particles: their bag's id; road segments: their road's
    const auto user_part = [&](int Family, std::int64_t P) {
        if (Family == kD3Airbag)
            return (P >= 0 && static_cast<std::size_t>(P) < bag_ids.size())
                       ? bag_ids[static_cast<std::size_t>(P)]
                       : P + 1;
        if (Family == kD3Road)
            return P + 1;
        return (P >= 0 && static_cast<std::size_t>(P) < part_ids.size())
                   ? part_ids[static_cast<std::size_t>(P)]
                   : P + 1;
    };
    std::vector<NDArray> eids, parts;
    // (kind: -1 a part, else the airbag or road family; id) -> (dim, cells)
    using D3Key = std::pair<int, std::int64_t>;
    std::map<D3Key, std::pair<int, std::vector<std::int64_t>>> by_part;
    std::vector<D3Key> part_order;
    std::int64_t base = 0;
    for (const D3Block& b : rCells.mBlocks) {
        const auto& ids = *family_ids[static_cast<std::size_t>(b.mFamily)];
        std::vector<std::int64_t> e(b.mElems.size()), p(b.mElems.size());
        const int dim = d3_dim(b.mType);
        const int kind = (b.mFamily == kD3Airbag || b.mFamily == kD3Road) ? b.mFamily : -1;
        for (std::size_t i = 0; i < b.mElems.size(); ++i) {
            e[i] = ids[b.mElems[i]];
            p[i] = user_part(b.mFamily, b.mParts[i]);
            const D3Key key{kind, p[i]};
            auto [it, inserted] = by_part.try_emplace(key, -1, std::vector<std::int64_t>{});
            if (inserted)
                part_order.push_back(key);
            it->second.first = std::max(it->second.first, dim);
            it->second.second.push_back(base + static_cast<std::int64_t>(i));
        }
        base += static_cast<std::int64_t>(b.mElems.size());
        eids.push_back(d3_int_array(e));
        parts.push_back(d3_int_array(p));
    }
    mesh.AddCellData("lsdyna:eid", std::move(eids));
    mesh.AddCellData("lsdyna:part", std::move(parts));
    for (const D3Key& key : part_order) {
        auto& [dim, members] = by_part[key];
        const std::int64_t pid = key.second;
        std::string name;
        if (key.first == kD3Airbag) {
            name = "Airbag " + std::to_string(pid);
        } else if (key.first == kD3Road) {
            name = "Rigid road " + std::to_string(pid);
        } else {
            const auto title = g.mPartTitles.find(pid);
            name = (title != g.mPartTitles.end() && !title->second.empty())
                       ? title->second
                       : "Part " + std::to_string(pid);
        }
        mesh.AddRegion(Region(name, RegionKind::Cell, dim, pid, d3_int_array(members)));
    }
    return mesh;
}

// --- one state -----------------------------------------------------------------------

// Named per-cell arrays assembled from the element families, NaN-padded to a
// common (points, width) across blocks.
class D3CellArrays {
public:
    explicit D3CellArrays(const D3Cells& rCells) : mrCells(rCells) {}

    // rRows: (elements, points * width) of Family, in element order.
    void Add(const std::string& rName, int Family, std::vector<double> Rows, std::size_t Points,
             std::size_t Width) {
        auto it = mEntries.find(rName);
        if (it == mEntries.end()) {
            it = mEntries.emplace(rName, Entry{}).first;
            mOrder.push_back(rName);
        }
        Entry& e = it->second;
        e.mPoints = std::max(e.mPoints, Points);
        e.mWidth = std::max(e.mWidth, Width);
        e.mParts.push_back({Family, std::move(Rows), Points, Width});
    }

    void Emit(Mesh& rMesh) const {
        for (const std::string& name : mOrder) {
            const Entry& e = mEntries.at(name);
            const std::size_t cols = e.mPoints * e.mWidth;
            std::vector<NDArray> out;
            for (const D3Block& b : mrCells.mBlocks) {
                NDArray a = cols == 1 ? NDArray(DType::Float64, {b.mElems.size()})
                                      : NDArray(DType::Float64, {b.mElems.size(), cols});
                std::fill(a.As<double>(), a.As<double>() + a.Size(), kD3Nan);
                out.push_back(std::move(a));
            }
            for (const Part& p : e.mParts) {
                const auto& where = mrCells.mWhere[static_cast<std::size_t>(p.mFamily)];
                for (std::size_t el = 0; el < where.size(); ++el) {
                    double* dst = out[where[el].first].As<double>() + where[el].second * cols;
                    const double* src = p.mRows.data() + el * p.mPoints * p.mWidth;
                    for (std::size_t q = 0; q < p.mPoints; ++q)
                        for (std::size_t c = 0; c < p.mWidth; ++c)
                            dst[q * e.mWidth + c] = src[q * p.mWidth + c];
                }
            }
            rMesh.AddCellData(name, std::move(out));
            if (e.mPoints > 1) {
                NDArray layout(DType::Int64, {2});
                layout.As<std::int64_t>()[0] = static_cast<std::int64_t>(e.mPoints);
                layout.As<std::int64_t>()[1] = static_cast<std::int64_t>(e.mWidth);
                rMesh.AddFieldData("lsdyna_d3plot:layout:" + name, std::move(layout));
            }
        }
    }

private:
    struct Part {
        int mFamily;
        std::vector<double> mRows;
        std::size_t mPoints, mWidth;
    };
    struct Entry {
        std::size_t mPoints = 0, mWidth = 0;
        std::vector<Part> mParts;
    };
    const D3Cells& mrCells;
    std::map<std::string, Entry> mEntries;
    std::vector<std::string> mOrder;
};

// Columns [C0, C0 + Count) of each row of a (rows, Stride) table, as (rows, Count).
std::vector<double> d3_columns(const double* pData, std::size_t Rows, std::size_t Stride,
                               std::size_t C0, std::size_t Count) {
    // Never past a row: columns beyond the stride are NaN.
    std::vector<double> out(Rows * Count, std::numeric_limits<double>::quiet_NaN());
    const std::size_t take = C0 < Stride ? std::min(Count, Stride - C0) : 0;
    for (std::size_t r = 0; r < Rows; ++r)
        std::copy(pData + r * Stride + C0, pData + r * Stride + C0 + take, out.data() + r * Count);
    return out;
}

// Per layered row (rows, Layers, LayerWidth) the columns [C0, C0 + Count) of
// every layer, as (rows, Layers * Count).
std::vector<double> d3_layer_columns(const double* pData, std::size_t Rows, std::size_t Stride,
                                     std::size_t Layers, std::size_t LayerWidth, std::size_t C0,
                                     std::size_t Count) {
    // The last layer's columns must end inside the row; otherwise the last
    // row's read runs past the state record.
    if (Layers > 0 && (Layers - 1) * LayerWidth + C0 + Count > Stride)
        d3_fail("a layered element record is wider than its row");
    std::vector<double> out(Rows * Layers * Count);
    for (std::size_t r = 0; r < Rows; ++r)
        for (std::size_t l = 0; l < Layers; ++l)
            std::copy(pData + r * Stride + l * LayerWidth + C0,
                      pData + r * Stride + l * LayerWidth + C0 + Count,
                      out.data() + (r * Layers + l) * Count);
    return out;
}

void d3_read_state(const D3File& rF, Mesh& rMesh, const D3Cells& rCells, std::size_t Index,
                   const ReadOptions& rOpts) {
    const D3Header& h = rF.mHeader;
    const D3Geometry& g = rF.mGeometry;
    const std::vector<double> s = rF.State(Index);
    const auto nn = static_cast<std::size_t>(h.mNodes);
    const auto want = [&](const std::string& rName) { return rOpts.WantsArray(rName); };
    std::size_t k = 1;
    const auto nglbv = static_cast<std::size_t>(std::max<std::int64_t>(h.mGlobals, 0));
    const double* globals = s.data() + k;
    k += nglbv;

    std::size_t gi = 0;
    for (const auto& [name, width] :
         std::array<std::pair<const char*, std::size_t>, 4>{{{"global_kinetic_energy", 1},
                                                             {"global_internal_energy", 1},
                                                             {"global_total_energy", 1},
                                                             {"global_velocity", 3}}}) {
        if (gi + width <= nglbv) {
            if (want(name)) {
                NDArray a(DType::Float64, {width});
                std::copy(globals + gi, globals + gi + width, a.As<double>());
                rMesh.AddFieldData(name, std::move(a));
            }
            gi += width;
        }
    }
    const auto nparts = static_cast<std::size_t>(std::max<std::int64_t>(
        h.R("nummat8") + h.R("nummat2") + h.R("nummat4") + h.R("nummatt") + g.mRigidBodies, 0));
    for (const auto& [name, width] :
         std::array<std::pair<const char*, std::size_t>, 5>{{{"part_internal_energy", 1},
                                                             {"part_kinetic_energy", 1},
                                                             {"part_velocity", 3},
                                                             {"part_mass", 1},
                                                             {"part_hourglass_energy", 1}}}) {
        if (gi + width * nparts <= nglbv) {
            if (want(name)) {
                NDArray a = width == 3 ? NDArray(DType::Float64, {nparts, 3})
                                       : NDArray(DType::Float64, {nparts});
                std::copy(globals + gi, globals + gi + width * nparts, a.As<double>());
                rMesh.AddFieldData(name, std::move(a));
            }
            gi += width * nparts;
        }
    }

    // Node arrays cover every point: NaN at the airbag particles and road nodes.
    const std::size_t npts = rMesh.NumPoints();
    for (const auto& [name, comps] : d3_node_vars(h)) {
        const double* v = s.data() + k;
        k += comps * nn;
        if (name == "coordinates") {
            if (want("displacement")) {
                NDArray a(DType::Float64, {npts, 3});
                double* out = a.As<double>();
                std::fill(out, out + a.Size(), kD3Nan);
                for (std::size_t i = 0; i < 3 * nn; ++i)
                    out[i] = v[i] - g.mCoords[i];
                rMesh.AddPointData("displacement", std::move(a));
            }
        } else if (want(name)) {
            NDArray a = comps > 1 ? NDArray(DType::Float64, {npts, comps})
                                  : NDArray(DType::Float64, {npts});
            std::fill(a.As<double>(), a.As<double>() + a.Size(), kD3Nan);
            std::copy(v, v + comps * nn, a.As<double>());
            rMesh.AddPointData(name, std::move(a));
        }
    }

    D3CellArrays cells(rCells);
    const auto put = [&](const std::string& rName, int Family, std::vector<double> Rows,
                         std::size_t Points, std::size_t Width) {
        if (want(rName))
            cells.Add(rName, Family, std::move(Rows), Points, Width);
    };
    const auto zu = [](std::int64_t N) {
        return static_cast<std::size_t>(std::max<std::int64_t>(N, 0));
    };

    if (h.mNt3d > 0) {
        const std::size_t n = zu(h.mSolids), nt = zu(h.mNt3d);
        put("thermal_variables", kD3Solid, std::vector<double>(s.data() + k, s.data() + k + n * nt),
            1, nt);
        k += n * nt;
    }

    if (h.mSolids > 0 && h.mNv3d > 0) {
        const std::size_t n = zu(h.mSolids), nv = zu(h.mNv3d);
        const std::size_t layers = zu(h.SolidLayers());
        const std::size_t lw = nv / layers;
        const double* d = s.data() + k;
        std::size_t i = 0;
        if (h.mSolidStress) {
            put("stress", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, i, 6), layers, 6);
            i += 6;
        }
        if (h.mSolidPstrain) {
            put("effective_plastic_strain", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, i, 1),
                layers, 1);
            i += 1;
        }
        // history = [i, i + neiph) of each layer, less what is split off it
        std::size_t h0 = i;
        std::size_t hn = std::min(zu(h.mNeiph), lw > i ? lw - i : 0);
        const std::size_t nstrain = h.mElementStrain ? 6 : 0;
        if (nstrain && hn >= nstrain) {
            put("strain", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, h0 + hn - nstrain, 6),
                layers, 6);
            hn -= nstrain;
        }
        if (h.mPlasticStrainTensor && hn >= 6) {
            put("plastic_strain_tensor", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, h0, 6),
                layers, 6);
            h0 += 6;
            hn -= 6;
        }
        if (h.mThermalStrainTensor && hn >= 6) {
            put("thermal_strain_tensor", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, h0, 6),
                layers, 6);
            h0 += 6;
            hn -= 6;
        }
        if (hn)
            put("history_variables", kD3Solid, d3_layer_columns(d, n, nv, layers, lw, h0, hn),
                layers, hn);
        k += n * nv;
    }

    if (h.mTshells > 0 && h.mNv3dt > 0) {
        const std::size_t n = zu(h.mTshells), nv = zu(h.mNv3dt), nl = zu(h.mLayers),
                          nh = zu(h.mNeips);
        const std::size_t lw = 6 * h.mShellStress + h.mShellPstrain + nh;
        const std::size_t nlayer = nl * lw;
        const double* d = s.data() + k;
        std::size_t i = 0;
        if (h.mShellStress) {
            put("stress", kD3Tshell, d3_layer_columns(d, n, nv, nl, lw, i, 6), nl, 6);
            i += 6;
        }
        if (h.mShellPstrain) {
            put("effective_plastic_strain", kD3Tshell, d3_layer_columns(d, n, nv, nl, lw, i, 1), nl,
                1);
            i += 1;
        }
        if (nh)
            put("history_variables", kD3Tshell, d3_layer_columns(d, n, nv, nl, lw, i, nh), nl, nh);
        if (h.mElementStrain) {
            put("strain_inner", kD3Tshell, d3_columns(d, n, nv, nlayer, 6), 1, 6);
            put("strain_outer", kD3Tshell, d3_columns(d, n, nv, nlayer + 6, 6), 1, 6);
        }
        k += n * nv;
    }

    if (h.mBeams > 0 && h.mNv1d > 0) {
        if (h.mNeipb < 0)
            d3_fail("negative beam history count");
        const std::size_t n = zu(h.mBeams), nv = zu(h.mNv1d), nh = zu(h.mNeipb);
        const std::int64_t nl_signed = (-3 * h.mNeipb + h.mNv1d - 6) / (h.mNeipb + 5);
        const std::size_t nl = zu(nl_signed);
        const double* d = s.data() + k;
        put("beam_axial_force", kD3Beam, d3_columns(d, n, nv, 0, 1), 1, 1);
        put("beam_shear_force", kD3Beam, d3_columns(d, n, nv, 1, 2), 1, 2);
        put("beam_bending_moment", kD3Beam, d3_columns(d, n, nv, 3, 2), 1, 2);
        put("beam_torsion_moment", kD3Beam, d3_columns(d, n, nv, 5, 1), 1, 1);
        if (nl > 0) {
            if (6 + 5 * nl > nv)
                d3_fail("a beam record is narrower than its integration points");
            const double* layered = d + 6;
            put("beam_axial_stress", kD3Beam, d3_layer_columns(layered, n, nv, nl, 5, 0, 1), nl, 1);
            put("beam_shear_stress", kD3Beam, d3_layer_columns(layered, n, nv, nl, 5, 1, 2), nl, 2);
            put("effective_plastic_strain", kD3Beam, d3_layer_columns(layered, n, nv, nl, 5, 3, 1),
                nl, 1);
            put("beam_axial_strain", kD3Beam, d3_layer_columns(layered, n, nv, nl, 5, 4, 1), nl, 1);
        }
        if (nh && nv > 6 + 5 * nl)
            put("history_variables", kD3Beam, d3_columns(d, n, nv, 6 + 5 * nl, nv - 6 - 5 * nl),
                3 + nl, nh);
        k += n * nv;
    }

    const std::int64_t nreduced = h.mShells - g.mRigidShells;
    if (nreduced > 0 && h.mNv2d > 0) {
        const std::size_t nv = zu(h.mNv2d), nl = zu(h.mLayers), nh = zu(h.mNeips);
        const std::size_t lw = 6 * h.mShellStress + h.mShellPstrain + nh;
        const std::size_t nlayer = nl * lw;
        const std::size_t n = zu(h.mShells);
        std::vector<double> full;
        const double* d = s.data() + k;
        k += zu(nreduced) * nv;
        if (g.mRigidShells) {
            full.assign(n * nv, kD3Nan);
            std::size_t r = 0;
            for (std::size_t e = 0; e < n; ++e) {
                const std::int64_t part = d3_zero_based(g.mShells[5 * e + 4]);
                const bool rigid = part >= 0 &&
                                   static_cast<std::size_t>(part) < g.mPartMattype.size() &&
                                   g.mPartMattype[static_cast<std::size_t>(part)] == 20;
                if (!rigid) {
                    std::copy(d + r * nv, d + (r + 1) * nv, full.data() + e * nv);
                    ++r;
                }
            }
            d = full.data();
        }
        std::size_t i = 0;
        if (h.mShellStress) {
            put("stress", kD3Shell, d3_layer_columns(d, n, nv, nl, lw, i, 6), nl, 6);
            i += 6;
        }
        if (h.mShellPstrain) {
            put("effective_plastic_strain", kD3Shell, d3_layer_columns(d, n, nv, nl, lw, i, 1), nl,
                1);
            i += 1;
        }
        if (nh)
            put("history_variables", kD3Shell, d3_layer_columns(d, n, nv, nl, lw, i, nh), nl, nh);
        std::size_t j = nlayer;
        if (h.mShellForces) {
            put("shell_bending_moment", kD3Shell, d3_columns(d, n, nv, j, 3), 1, 3);
            put("shell_shear_force", kD3Shell, d3_columns(d, n, nv, j + 3, 2), 1, 2);
            put("shell_normal_force", kD3Shell, d3_columns(d, n, nv, j + 5, 3), 1, 3);
            j += 8;
        }
        if (h.mShellExtra) {
            put("thickness", kD3Shell, d3_columns(d, n, nv, j, 1), 1, 1);
            put("shell_element_variables", kD3Shell, d3_columns(d, n, nv, j + 1, 2), 1, 2);
            j += 3;
        }
        if (h.mElementStrain) {
            put("strain_inner", kD3Shell, d3_columns(d, n, nv, j, 6), 1, 6);
            put("strain_outer", kD3Shell, d3_columns(d, n, nv, j + 6, 6), 1, 6);
            j += 12;
        }
        if (h.mShellExtra) {
            put("internal_energy", kD3Shell, d3_columns(d, n, nv, j, 1), 1, 1);
            j += 1;
        }
        if (h.mPlasticStrainTensor) {
            // Per layer; a file can hold it at fewer points than it has layers
            // (a composite shell's 10 layers, 3 tensors): as many as fit.
            const std::size_t points = std::min(nl, j < nv ? (nv - j) / 6 : 0);
            if (points < nl)
                log::warn("LS-DYNA d3plot: the shells' plastic strain tensor has {} of {} layers",
                          points, nl);
            if (points > 0)
                put("plastic_strain_tensor", kD3Shell, d3_columns(d, n, nv, j, 6 * points), points,
                    6);
            j += 6 * points;
        }
        if (h.mThermalStrainTensor && j + 6 <= nv) {
            put("thermal_strain_tensor", kD3Shell, d3_columns(d, n, nv, j, 6), 1, 6);
            j += 6;
        }
    }

    // deletion (before the SPH data)
    std::array<std::size_t, kD3Families> alive_at{};
    std::array<bool, kD3Families> has_alive{};
    if (h.mNodeDeletion) {
        if (want("lsdyna:alive")) {
            NDArray a(DType::Int8, {npts});
            std::fill(a.As<std::int8_t>(), a.As<std::int8_t>() + npts, std::int8_t{0});
            for (std::size_t i = 0; i < nn; ++i)
                a.As<std::int8_t>()[i] = s[k + i] != 0.0 ? 1 : 0;
            rMesh.AddPointData("lsdyna:alive", std::move(a));
        }
        k += nn;
    } else if (h.mElementDeletion) {
        for (const auto& [family, count] :
             std::array<std::pair<int, std::int64_t>, 4>{{{kD3Solid, h.mSolids},
                                                          {kD3Tshell, h.mTshells},
                                                          {kD3Shell, h.mShells},
                                                          {kD3Beam, h.mBeams}}}) {
            alive_at[static_cast<std::size_t>(family)] = k;
            has_alive[static_cast<std::size_t>(family)] = true;
            k += zu(count);
        }
    }

    // SPH particles: a material word (negative once deleted), then the
    // variables ISPHFG(2..11) flag, each as many words as its flag says.
    if (h.mSph > 0 && g.mSphVars > 0) {
        const std::size_t n = zu(h.mSph), nv = zu(g.mSphVars);
        const double* d = s.data() + k;
        alive_at[kD3Sph] = k;
        has_alive[kD3Sph] = true;
        const auto& f = g.mSphFlags;
        const std::size_t strain = static_cast<std::size_t>(std::abs(f[8]));
        const std::size_t strain6 = std::min<std::size_t>(strain, 6);
        const std::array<std::pair<const char*, std::size_t>, 11> vars{{
            {"sph_radius", zu(f[1])},
            {"sph_pressure", zu(f[2])},
            {"stress", zu(f[3])},
            {"effective_plastic_strain", zu(f[4])},
            {"density", zu(f[5])},
            {"internal_energy", zu(f[6])},
            {"sph_neighbors", zu(f[7])},
            {"strain", strain6},
            {"strain_rate", strain - strain6},
            {"mass", zu(f[9])},
            {"history_variables", zu(f[10])},
        }};
        std::size_t i = 1;
        for (const auto& [name, width] : vars) {
            if (width > 0)
                put(name, kD3Sph, d3_columns(d, n, nv, i, width), 1, width);
            i += width;
        }
        k += n * nv;
    }

    bool any_alive = false;
    for (bool b : has_alive)
        any_alive = any_alive || b;
    if (want("lsdyna:alive") && !rCells.mBlocks.empty() && any_alive) {
        std::vector<NDArray> out;
        for (const D3Block& b : rCells.mBlocks) {
            NDArray a(DType::Int8, {b.mElems.size()});
            const auto fam = static_cast<std::size_t>(b.mFamily);
            for (std::size_t i = 0; i < b.mElems.size(); ++i) {
                std::int8_t v = 1;  // families without deletion data live
                if (has_alive[fam]) {
                    const std::size_t at =
                        alive_at[fam] + b.mElems[i] * (fam == kD3Sph ? zu(g.mSphVars) : 1);
                    v = fam == kD3Sph ? (s[at] >= 0.0 ? 1 : 0) : (s[at] != 0.0 ? 1 : 0);
                }
                a.As<std::int8_t>()[i] = v;
            }
            out.push_back(std::move(a));
        }
        rMesh.AddCellData("lsdyna:alive", std::move(out));
    }

    // airbags: per bag its state variables, then per particle its variables,
    // named in the geometry section (positions become the particles' points)
    if (g.mHasAirbag) {
        const std::size_t ngeom = zu(g.mAirbagGeom), nvar = zu(g.mAirbagVar);
        const std::size_t npart = zu(g.mAirbagParticles), nst = zu(g.mAirbagStateGeom);
        const std::size_t nbags = zu(h.mAirbags);
        const std::vector<std::int64_t> ints = rF.StateInts(Index);
        const std::size_t bag0 = k, part0 = k + nbags * nst;
        k = part0 + npart * nvar;
        const auto value = [&](std::size_t At, std::int64_t Type) {
            // type 1: an integer stored in the word
            return Type == 1 ? static_cast<double>(ints[At]) : s[At];
        };
        const auto airbag_name = [](std::string Name) {
            // lower case, words joined by underscores ("Bag Vol" -> "bag_vol")
            std::string out;
            bool gap = false;
            for (char c : Name) {
                if (c == ' ' || c == '\t') {
                    gap = !out.empty();
                    continue;
                }
                if (gap)
                    out.push_back('_');
                gap = false;
                out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
            }
            return out;
        };
        for (std::size_t j = 0; j < nst; ++j) {
            const std::string key =
                "lsdyna:airbag:" + airbag_name(g.mAirbagNames[ngeom + nvar + j]);
            if (!want(key))
                continue;
            NDArray a(DType::Float64, {nbags});
            for (std::size_t b = 0; b < nbags; ++b)
                a.As<double>()[b] = value(bag0 + b * nst + j, g.mAirbagTypes[ngeom + nvar + j]);
            rMesh.AddFieldData(key, std::move(a));
        }
        const auto column_of = [&](const char* pName) -> std::size_t {
            for (std::size_t j = 0; j < nvar; ++j)
                if (g.mAirbagNames[ngeom + j] == pName)
                    return j;
            return kD3Npos;
        };
        const std::size_t px = column_of("Pos x"), py = column_of("Pos y"), pz = column_of("Pos z");
        if (px != kD3Npos && py != kD3Npos && pz != kD3Npos) {
            NDArray pts = rMesh.Points();
            double* p = pts.As<double>();
            for (std::size_t i = 0; i < npart; ++i) {
                p[3 * (nn + i)] = s[part0 + i * nvar + px];
                p[3 * (nn + i) + 1] = s[part0 + i * nvar + py];
                p[3 * (nn + i) + 2] = s[part0 + i * nvar + pz];
            }
            rMesh.AssignPoints(std::move(pts));
        }
        const std::size_t vx = column_of("Vel x"), vy = column_of("Vel y"), vz = column_of("Vel z");
        if (vx != kD3Npos && vy != kD3Npos && vz != kD3Npos) {
            std::vector<double> rows(3 * npart);
            for (std::size_t i = 0; i < npart; ++i) {
                rows[3 * i] = s[part0 + i * nvar + vx];
                rows[3 * i + 1] = s[part0 + i * nvar + vy];
                rows[3 * i + 2] = s[part0 + i * nvar + vz];
            }
            put("velocity", kD3Airbag, std::move(rows), 1, 3);
        }
        for (std::size_t j = 0; j < nvar; ++j) {
            const std::string& name = g.mAirbagNames[ngeom + j];
            if (name.rfind("Pos ", 0) == 0 || name.rfind("Vel ", 0) == 0)
                continue;
            std::vector<double> rows(npart);
            for (std::size_t i = 0; i < npart; ++i)
                rows[i] = value(part0 + i * nvar + j, g.mAirbagTypes[ngeom + j]);
            put("airbag_" + airbag_name(name), kD3Airbag, std::move(rows), 1, 1);
        }
    }

    // rigid roads: per road its displacement and velocity
    if (g.mRoads > 0) {
        const std::size_t nr = zu(g.mRoads);
        for (const auto& [name, j] : std::array<std::pair<const char*, std::size_t>, 2>{
                 {{"lsdyna:road_displacement", 0}, {"lsdyna:road_velocity", 1}}}) {
            if (!want(name))
                continue;
            NDArray a(DType::Float64, {nr, 3});
            for (std::size_t r = 0; r < nr; ++r)
                for (std::size_t c = 0; c < 3; ++c)
                    a.As<double>()[3 * r + c] = s[k + 6 * r + 3 * j + c];
            rMesh.AddFieldData(name, std::move(a));
        }
        k += 6 * nr;
    }

    // rigid bodies: centre of mass, rotation matrix and (unless reduced)
    // velocity, rotational velocity, acceleration, rotational acceleration
    if (h.mRigidBodies && g.mRigidBodyMotions > 0) {
        const std::size_t nr = zu(g.mRigidBodyMotions), nv = h.mReducedRigidBodies ? 12 : 24;
        std::vector<std::pair<const char*, std::size_t>> fields{{"coordinates", 3},
                                                                {"rotation", 9}};
        if (!h.mReducedRigidBodies)
            fields.insert(fields.end(), {{"velocity", 3},
                                         {"rotational_velocity", 3},
                                         {"acceleration", 3},
                                         {"rotational_acceleration", 3}});
        std::size_t i = 0;
        for (const auto& [name, width] : fields) {
            const std::string key = std::string("lsdyna:rigid_body_") + name;
            if (want(key)) {
                NDArray a(DType::Float64, {nr, width});
                for (std::size_t r = 0; r < nr; ++r)
                    std::copy(s.data() + k + r * nv + i, s.data() + k + r * nv + i + width,
                              a.As<double>() + r * width);
                rMesh.AddFieldData(key, std::move(a));
            }
            i += width;
        }
        k += nr * nv;
    }
    cells.Emit(rMesh);
}

}  // namespace

bool is_d3plot_filename(const std::string& rPath) {
    std::string name = fs::path(rPath).filename().string();
    for (char& c : name)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return name == "d3plot" || name == "d3part";
}

bool is_d3plot_head(const char* pHead, std::size_t Size) {
    const auto sniffed = d3_sniff(pHead, Size);
    if (!sniffed)
        return false;
    const D3Words w{pHead, Size, sniffed->first, sniffed->second};
    const double version = w.FloatAt(14);
    if (!(version == 0.0 || (version >= 900.0 && version < 100000.0)))
        return false;
    for (int i : {16, 18, 28, 31, 40})
        if (w.IntAt(static_cast<std::size_t>(i)) < 0)
            return false;
    for (std::size_t i = 0; i < 40; ++i) {
        const auto b = static_cast<unsigned char>(pHead[i]);
        if (b != 0 && (b < 32 || b >= 127))
            return false;
    }
    return true;
}

Mesh read_lsdyna_d3plot(const std::string& rPath, const ReadOptions& rOpts) {
    const D3File file(rPath);
    D3Cells cells;
    Mesh mesh = d3_build_mesh(file, cells);
    const std::size_t n = file.mStates.size();
    if (n == 0) {
        if (rOpts.mTimeStep != 0 && rOpts.mTimeStep != -1)
            throw ReadError("time step " + std::to_string(rOpts.mTimeStep) +
                            " is out of range: the file has no states");
        return mesh;
    }
    const std::size_t index = rOpts.ResolveTimeStep(n);
    const std::vector<double> times = file.Times();
    mesh.AddFieldData(kSequenceTimeKey, d3_scalar(times[index], DType::Float64));
    mesh.AddFieldData("lsdyna:state", d3_scalar(static_cast<double>(index), DType::Int64));
    if (rOpts.mPointsOnly)
        return mesh;
    d3_read_state(file, mesh, cells, index, rOpts);
    return mesh;
}

std::vector<double> lsdyna_d3plot_time_values(const std::string& rPath) {
    return D3File(rPath).Times();
}

MeshMetadata read_lsdyna_d3plot_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    ReadOptions options = rOpts;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(read_lsdyna_d3plot(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "lsdyna_d3plot";
    meta.mTimeValues = lsdyna_d3plot_time_values(rPath);
    return meta;
}

// --- writing ------------------------------------------------------------------------------

namespace {

constexpr std::size_t kD3wDefaultFileBytes = std::size_t{1} << 30;
constexpr double kD3wVersion = 971.0;
enum D3wFamily : int { kD3wSolid = 0, kD3wBeam = 1, kD3wShell = 2, kD3wSph = 3 };

[[noreturn]] void d3w_fail(const std::string& rMessage) {
    throw WriteError("LS-DYNA d3plot: " + rMessage);
}

/// Little-endian words of `Ws` bytes (4, or 8 in double precision).
void d3w_i32(std::string& rOut, std::int64_t V, int Ws = 4) {
    if (Ws == 8) {
        const auto u = static_cast<std::uint64_t>(V);
        for (int k = 0; k < 8; ++k)
            rOut.push_back(static_cast<char>((u >> (8 * k)) & 0xff));
        return;
    }
    const auto u = static_cast<std::uint32_t>(static_cast<std::int32_t>(V));
    for (int k = 0; k < 4; ++k)
        rOut.push_back(static_cast<char>((u >> (8 * k)) & 0xff));
}

void d3w_f32(std::string& rOut, double V, int Ws = 4) {
    if (Ws == 8) {
        const auto u = std::bit_cast<std::uint64_t>(V);
        for (int k = 0; k < 8; ++k)
            rOut.push_back(static_cast<char>((u >> (8 * k)) & 0xff));
        return;
    }
    const auto u = std::bit_cast<std::uint32_t>(static_cast<float>(V));
    for (int k = 0; k < 4; ++k)
        rOut.push_back(static_cast<char>((u >> (8 * k)) & 0xff));
}

int d3w_family(std::string_view Type) {
    if (Type == "hexahedron" || Type == "wedge" || Type == "pyramid" || Type == "tetra" ||
        Type == "tetra10" || Type == "hexahedron20" || Type == "hexahedron27")
        return kD3wSolid;
    if (Type == "line")
        return kD3wBeam;
    if (Type == "triangle" || Type == "quad" || Type == "quad8")
        return kD3wShell;
    if (Type == "vertex")
        return kD3wSph;
    return -1;
}

double d3w_time_of(const Mesh& rMesh) {
    return rMesh.HasFieldData("meshio:time") && rMesh.FieldData("meshio:time").Size() > 0
               ? detail::read_double(rMesh.FieldData("meshio:time"), 0)
               : 0.0;
}

/// The default components of an element variable (0: from the layout array).
int d3w_default_comps(const std::string& rName) {
    static const std::unordered_map<std::string, int> table = {
        {"stress", 6},
        {"effective_plastic_strain", 1},
        {"plastic_strain_tensor", 6},
        {"thermal_strain_tensor", 6},
        {"strain", 6},
        {"thermal_variables", 0},
        {"history_variables", 0},
        {"shell_bending_moment", 3},
        {"shell_shear_force", 2},
        {"shell_normal_force", 3},
        {"thickness", 1},
        {"shell_element_variables", 2},
        {"strain_inner", 6},
        {"strain_outer", 6},
        {"internal_energy", 1},
        {"beam_axial_force", 1},
        {"beam_shear_force", 2},
        {"beam_bending_moment", 2},
        {"beam_torsion_moment", 1},
        {"beam_axial_stress", 1},
        {"beam_shear_stress", 2},
        {"beam_axial_strain", 1},
        {"sph_radius", 1},
        {"sph_pressure", 1},
        {"density", 1},
        {"sph_neighbors", 1},
        {"strain_rate", 6},
        {"mass", 1},
    };
    const auto it = table.find(rName);
    return it == table.end() ? 0 : it->second;
}

/// One block's values of a cell array as (cells, points, comps); mPoints 0 when unusable.
struct D3wBlockArray {
    std::size_t mPoints = 0;
    std::size_t mComps = 0;
    std::vector<double> mV;
};

std::vector<D3wBlockArray> d3w_cell_array(const Mesh& rMesh, const std::string& rName,
                                          bool& rPresent) {
    rPresent = rMesh.HasCellData(rName);
    std::vector<D3wBlockArray> out(rMesh.NumCellBlocks());
    if (!rPresent)
        return out;
    std::int64_t layout0 = 0, layout1 = 0;
    bool has_layout = false;
    const std::string key = "lsdyna_d3plot:layout:" + rName;
    if (rMesh.HasFieldData(key) && rMesh.FieldData(key).Size() >= 2) {
        layout0 = detail::read_int(rMesh.FieldData(key), 0);
        layout1 = detail::read_int(rMesh.FieldData(key), 1);
        has_layout = true;
    }
    const auto dflt = static_cast<std::size_t>(d3w_default_comps(rName));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::size_t n = rMesh.Cells(b).NumCells();
        if (b >= rMesh.CellDataNumBlocks(rName))
            continue;
        const NDArray& a = rMesh.CellData(rName, b);
        if (n == 0 || a.Size() % n != 0)
            continue;
        const std::size_t width = a.Size() / n;
        std::size_t points = 0, comps = 0;
        if (has_layout && layout0 > 0 && layout1 > 0 &&
            width == static_cast<std::size_t>(layout0) * static_cast<std::size_t>(layout1)) {
            points = static_cast<std::size_t>(layout0);
            comps = static_cast<std::size_t>(layout1);
        } else if (dflt && width > 0 && width % dflt == 0) {
            points = width / dflt;
            comps = dflt;
        } else if (!dflt && width > 0) {
            points = 1;
            comps = width;
        } else {
            continue;
        }
        out[b].mPoints = points;
        out[b].mComps = comps;
        out[b].mV.resize(a.Size());
        for (std::size_t i = 0; i < a.Size(); ++i)
            out[b].mV[i] = detail::read_double(a, i);
    }
    return out;
}

struct D3wRange {
    std::size_t mBlock, mFirst, mCount;
};

/// The geometry of a file and the layout of its states.
class D3wModel {
public:
    explicit D3wModel(const Mesh& rMesh) {
        mN = rMesh.NumPoints();
        if (mN == 0)
            d3w_fail("the mesh has no points");
        const NDArray& pts = rMesh.Points();
        const std::size_t width = pts.Size() / mN;
        mPoints.assign(mN * 3, 0.0);
        for (std::size_t i = 0; i < mN; ++i)
            for (std::size_t k = 0; k < 3 && k < width; ++k)
                mPoints[i * 3 + k] = detail::read_double(pts, i * width + k);
        mNid.resize(mN);
        if (rMesh.HasPointData("lsdyna:nid") && rMesh.PointData("lsdyna:nid").Size() == mN)
            for (std::size_t i = 0; i < mN; ++i)
                mNid[i] = detail::read_int(rMesh.PointData("lsdyna:nid"), i);
        else
            for (std::size_t i = 0; i < mN; ++i)
                mNid[i] = static_cast<std::int64_t>(i) + 1;

        std::size_t ncells = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            ncells += rMesh.Cells(b).NumCells();
        const std::vector<std::int64_t> part_of_cell = PartAssignment(rMesh, ncells);

        std::map<std::string, std::size_t> skipped;
        std::size_t airbag = 0;
        std::size_t base = 0;
        const bool has_eid = rMesh.HasCellData("lsdyna:eid");
        const bool has_part = rMesh.HasCellData("lsdyna:part");
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const std::string type(cb.Type());
            const std::size_t count = cb.NumCells();
            const std::size_t k = cb.NodesPerCell();
            const int family = d3w_family(type);
            if (family < 0) {
                if (count)
                    skipped[type] += count;
                base += count;
                continue;
            }
            const detail::Int64View conn(cb.Conn());
            for (std::size_t i = 0; i < count * k; ++i)
                if (conn.Data()[i] < 0 || static_cast<std::size_t>(conn.Data()[i]) >= mN)
                    d3w_fail("a '" + type + "' cell refers to a point the mesh does not have");
            if (family == kD3wSph) {
                // an airbag particle reads back on a point with no coordinates yet
                bool airbag_block = false;
                for (std::size_t i = 0; i < count && !airbag_block; ++i)
                    for (std::size_t c = 0; c < 3; ++c)
                        airbag_block =
                            airbag_block ||
                            std::isnan(mPoints[static_cast<std::size_t>(conn.Data()[i]) * 3 + c]);
                if (airbag_block) {
                    airbag += count;
                    base += count;
                    continue;
                }
            }
            const std::size_t first = mRows[family].size();
            for (std::size_t r = 0; r < count; ++r) {
                const std::int64_t* row = conn.Data() + r * k;
                std::vector<std::int64_t> nodes;
                if (family == kD3wSolid && type == "tetra10") {
                    nodes.assign(row, row + 8);
                    mExtra.push_back({mRows[family].size(), {row[8], row[9]}});
                } else if (family == kD3wSolid && type == "hexahedron20") {
                    nodes.assign(row, row + 8);
                    mExtra20.push_back(
                        {mRows[family].size(), std::vector<std::int64_t>(row + 8, row + 20)});
                } else if (family == kD3wSolid && type == "hexahedron27") {
                    nodes.assign(row, row + 8);
                    mExtra27.push_back(
                        {mRows[family].size(), std::vector<std::int64_t>(row + 8, row + 27)});
                } else if (family == kD3wSolid) {
                    const auto eight = detail::expand_brick(type, row);
                    nodes.assign(eight.begin(), eight.end());
                } else {
                    nodes.assign(row, row + k);
                    if (type == "triangle") {
                        nodes.push_back(nodes[2]);
                    } else if (type == "quad8") {
                        mExtra8.push_back(
                            {mRows[family].size(),
                             std::vector<std::int64_t>(nodes.begin() + 4, nodes.end())});
                        nodes.resize(4);
                    }
                }
                mRows[family].push_back(std::move(nodes));
                if (has_part && rMesh.CellData("lsdyna:part", b).Size() > r)
                    mPartsOf[family].push_back(
                        detail::read_int(rMesh.CellData("lsdyna:part", b), r));
                else
                    mPartsOf[family].push_back(part_of_cell[base + r]);
                if (has_eid && rMesh.CellData("lsdyna:eid", b).Size() > r)
                    mEid[family].push_back(detail::read_int(rMesh.CellData("lsdyna:eid", b), r));
                else
                    mEid[family].push_back(static_cast<std::int64_t>(mRows[family].size()));
            }
            mRanges[family].push_back({b, first, count});
            base += count;
        }
        for (const auto& [type, count] : skipped)
            log::warn("LS-DYNA d3plot: {} '{}' cell(s) are not written", count, type);
        if (airbag)
            log::warn("LS-DYNA d3plot: {} airbag particle(s) are not written", airbag);
        for (int f = 0; f < 4; ++f)
            mCounts[f] = mRows[f].size();
        if (!mCounts[0] && !mCounts[1] && !mCounts[2] && !mCounts[3])
            d3w_fail("the mesh has no solid, beam, shell or SPH cells");
        // rigid bodies: their parts (user ids) from the reader's field data
        if (rMesh.HasFieldData("lsdyna:rigid_body_part")) {
            const NDArray& rp = rMesh.FieldData("lsdyna:rigid_body_part");
            for (std::size_t i = 0; i < rp.Size(); ++i)
                mRigidParts.push_back(detail::read_int(rp, i));
        }
        mNrigid = mRigidParts.size();
        // every part, sorted; each family's number of distinct parts (SPH particles count with
        // the solids, as LS-DYNA counts them)
        std::set<std::int64_t> all(mRigidParts.begin(), mRigidParts.end());
        for (int f = 0; f < 4; ++f)
            all.insert(mPartsOf[f].begin(), mPartsOf[f].end());
        {
            std::set<std::int64_t> solid_sph(mPartsOf[kD3wSolid].begin(),
                                             mPartsOf[kD3wSolid].end());
            solid_sph.insert(mPartsOf[kD3wSph].begin(), mPartsOf[kD3wSph].end());
            mNummat[kD3wSolid] = solid_sph.size();
            mNummat[kD3wBeam] =
                std::set<std::int64_t>(mPartsOf[kD3wBeam].begin(), mPartsOf[kD3wBeam].end()).size();
            mNummat[kD3wShell] =
                std::set<std::int64_t>(mPartsOf[kD3wShell].begin(), mPartsOf[kD3wShell].end())
                    .size();
        }
        mNparts = mNummat[0] + mNummat[1] + mNummat[2] + mNrigid;
        mPartIds.assign(all.begin(), all.end());
        for (std::size_t i = 0; i < mPartIds.size(); ++i)
            mPartIndex[mPartIds[i]] = static_cast<std::int64_t>(i) + 1;
        mSignature = std::to_string(mN);
        for (const std::int64_t part : mRigidParts)
            mSignature += "r" + std::to_string(part);
        for (int f = 0; f < 4; ++f) {
            mSignature += "|" + std::to_string(f);
            for (const auto& row : mRows[f]) {
                mSignature += ";";
                for (const std::int64_t v : row)
                    mSignature += std::to_string(v) + ",";
            }
        }
        for (const auto& [e, nodes] : mExtra) {
            mSignature += "|t" + std::to_string(e);
            for (const std::int64_t v : nodes)
                mSignature += "," + std::to_string(v);
        }
        for (const auto* table : {&mExtra20, &mExtra27, &mExtra8})
            for (const auto& [e, nodes] : *table) {
                mSignature += "|x" + std::to_string(e);
                for (const std::int64_t v : nodes)
                    mSignature += "," + std::to_string(v);
            }
    }

    std::size_t mN = 0;
    std::vector<double> mPoints;
    std::vector<std::int64_t> mNid;
    std::array<std::vector<D3wRange>, 4> mRanges;
    std::array<std::vector<std::vector<std::int64_t>>, 4> mRows;
    std::vector<std::pair<std::size_t, std::array<std::int64_t, 2>>> mExtra;
    /// hexahedron20 / hexahedron27 / quad8: the element and its nodes past the corners.
    std::vector<std::pair<std::size_t, std::vector<std::int64_t>>> mExtra20, mExtra27, mExtra8;
    std::array<std::vector<std::int64_t>, 4> mPartsOf, mEid;
    std::array<std::size_t, 4> mCounts{};
    std::array<std::size_t, 3> mNummat{};
    std::vector<std::int64_t> mRigidParts;  ///< each rigid body's part (user id)
    std::size_t mNrigid = 0, mNparts = 0;
    std::vector<std::int64_t> mPartIds;
    std::map<std::int64_t, std::int64_t> mPartIndex;
    std::map<std::int64_t, std::string> mTitles;
    std::string mSignature;

private:
    std::vector<std::int64_t> PartAssignment(const Mesh& rMesh, std::size_t NumCells) {
        std::set<std::int64_t> taken;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const Region& r = rMesh.Region(i);
            if (r.mKind == RegionKind::Cell && r.mTag > 0)
                taken.insert(r.mTag);
        }
        std::int64_t fresh = 1;
        std::vector<std::int64_t> part(NumCells, -1);
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const Region& r = rMesh.Region(i);
            if (r.mKind != RegionKind::Cell)
                continue;
            std::int64_t pid = r.mTag;
            if (pid <= 0) {
                while (taken.count(fresh))
                    ++fresh;
                pid = fresh;
                taken.insert(pid);
            }
            mTitles.emplace(pid, r.mName);
            for (std::size_t q = 0; q < r.NumEntries(); ++q) {
                const std::int64_t e = r.Entries()[q];
                if (e >= 0 && static_cast<std::size_t>(e) < NumCells && part[e] < 0)
                    part[e] = pid;
            }
        }
        const std::int64_t rest = taken.empty() ? 1 : *taken.rbegin() + 1;
        for (auto& p : part)
            if (p < 0)
                p = rest;
        return part;
    }
};

/// A family's element-variable columns: the `mLayer` slots (name, components) repeat for each of
/// `mLayers` layers, the `mTail` slots (name, components, count) follow once.
struct D3wLayout {
    std::size_t mLayers = 1;
    std::vector<std::pair<std::string, std::size_t>> mLayer;
    struct Tail {
        std::string mName;
        std::size_t mComps, mCount;
    };
    std::vector<Tail> mTail;

    std::size_t LayerWidth() const {
        std::size_t w = 0;
        for (const auto& slot : mLayer)
            w += slot.second;
        return w;
    }
    std::size_t Width() const {
        std::size_t w = mLayers * LayerWidth();
        for (const Tail& t : mTail)
            w += t.mComps * t.mCount;
        return w;
    }
};

/// What the arrays of all steps hold for one family: name -> {points, comps}.
using D3wSeen = std::map<std::string, std::array<std::size_t, 2>>;

const std::vector<std::pair<const char*, std::size_t>>& d3w_globals() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"global_kinetic_energy", 1},
        {"global_internal_energy", 1},
        {"global_total_energy", 1},
        {"global_velocity", 3}};
    return v;
}

const std::vector<std::pair<const char*, std::size_t>>& d3w_part_globals() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"part_internal_energy", 1},
        {"part_kinetic_energy", 1},
        {"part_velocity", 3},
        {"part_mass", 1},
        {"part_hourglass_energy", 1}};
    return v;
}

const std::vector<std::pair<const char*, std::size_t>>& d3w_beam_core() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {{"beam_axial_force", 1},
                                                                       {"beam_shear_force", 2},
                                                                       {"beam_bending_moment", 2},
                                                                       {"beam_torsion_moment", 1}};
    return v;
}

/// Per beam layer: axial stress, shear stress (2), plastic strain, axial strain.
const std::vector<std::pair<const char*, std::size_t>>& d3w_beam_layer() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"beam_axial_stress", 1},
        {"beam_shear_stress", 2},
        {"effective_plastic_strain", 1},
        {"beam_axial_strain", 1}};
    return v;
}

const std::vector<std::pair<const char*, std::size_t>>& d3w_shell_forces() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"shell_bending_moment", 3}, {"shell_shear_force", 2}, {"shell_normal_force", 3}};
    return v;
}

const std::vector<std::pair<const char*, std::size_t>>& d3w_shell_extra() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"thickness", 1}, {"shell_element_variables", 2}};
    return v;
}

const std::vector<std::pair<const char*, std::size_t>>& d3w_rigid_fields() {
    static const std::vector<std::pair<const char*, std::size_t>> v = {
        {"coordinates", 3},         {"rotation", 9},     {"velocity", 3},
        {"rotational_velocity", 3}, {"acceleration", 3}, {"rotational_acceleration", 3}};
    return v;
}

/// SPH particle variables in state order: name, flag index, components.
struct D3wSphVar {
    const char* mName;
    std::size_t mFlag, mComps;
};
const std::vector<D3wSphVar>& d3w_sph_vars() {
    static const std::vector<D3wSphVar> v = {
        {"sph_radius", 1, 1},   {"sph_pressure", 2, 1},
        {"stress", 3, 6},       {"effective_plastic_strain", 4, 1},
        {"density", 5, 1},      {"internal_energy", 6, 1},
        {"sph_neighbors", 7, 1}};
    return v;
}

/// The header: what the states of all steps need.
class D3wPlan {
public:
    D3wPlan(const D3wModel& rModel, const std::vector<const Mesh*>& rMeshes, bool Double)
        : mModel(rModel), mWs(Double ? 8 : 4) {
        const std::size_t nparts = rModel.mNparts;
        static const char* kNodal[] = {
            "displacement",         "velocity",        "acceleration",
            "temperature",          "heat_flux",       "mass_scaling",
            "temperature_gradient", "residual_forces", "residual_moments"};
        std::array<D3wSeen, 4> seen;
        for (const Mesh* pMesh : rMeshes) {
            const Mesh& mesh = *pMesh;
            for (const char* name : kNodal)
                if (mesh.HasPointData(name)) {
                    const std::size_t comps =
                        mesh.PointData(name).Size() / std::max<std::size_t>(rModel.mN, 1);
                    mNode[name] = std::max(mNode[name], comps);
                }
            mAlive = mAlive || mesh.HasCellData("lsdyna:alive");
            mAliveNodes = mAliveNodes || mesh.HasPointData("lsdyna:alive");
            for (const auto& g : d3w_globals())
                mGlobals = mGlobals || mesh.HasFieldData(g.first);
            // the per-part arrays may cover more parts than the mesh uses
            const std::size_t count = mesh.HasFieldData(d3w_part_globals()[0].first)
                                          ? mesh.FieldData(d3w_part_globals()[0].first).Size()
                                          : 0;
            if (count >= nparts && !mPartGlobals) {
                bool all = true;
                for (const auto& [name, width] : d3w_part_globals())
                    all = all && mesh.HasFieldData(name) &&
                          mesh.FieldData(name).Size() == width * count;
                if (all) {
                    mPartGlobals = true;
                    mPartCount = count;
                }
            }
            if (rModel.mNrigid)
                for (const auto& [name, width] : d3w_rigid_fields()) {
                    const std::string key = std::string("lsdyna:rigid_body_") + name;
                    if (mesh.HasFieldData(key) &&
                        mesh.FieldData(key).Size() == rModel.mNrigid * width)
                        mRigid[name] = width;
                }
            for (const char* name : kElementNames()) {
                bool present = false;
                const auto arrays = d3w_cell_array(mesh, name, present);
                if (!present)
                    continue;
                for (const int family : {kD3wSolid, kD3wBeam, kD3wShell, kD3wSph})
                    for (const D3wRange& range : rModel.mRanges[family]) {
                        const D3wBlockArray& a = arrays[range.mBlock];
                        if (a.mPoints == 0)
                            continue;
                        if (family == kD3wSolid && std::string(name) != "thermal_variables" &&
                            a.mPoints != 1 && a.mPoints != 8)
                            continue;
                        auto& entry = seen[family][name];
                        entry[0] = std::max(entry[0], a.mPoints);
                        entry[1] = std::max(entry[1], a.mComps);
                    }
            }
        }
        mSeen = std::move(seen);
        mGlobals = mGlobals || mPartGlobals;
        // the header's material counts add up to the parts the arrays cover
        mNummat = rModel.mNummat;
        if (mPartGlobals)
            mNummat[kD3wSolid] += mPartCount - nparts;
        mNglbv = mGlobals ? 6 + (mPartGlobals ? 7 * mPartCount : 0) : 0;
        ElementLayouts();
        NodalFlags();
        mWords = StateWords();
        if (mWords < 2) {
            mDisp = true;
            mWords = StateWords();
        }
    }

    static const std::vector<const char*>& kElementNames() {
        static const std::vector<const char*> v = {"stress",
                                                   "effective_plastic_strain",
                                                   "plastic_strain_tensor",
                                                   "thermal_strain_tensor",
                                                   "strain",
                                                   "thermal_variables",
                                                   "history_variables",
                                                   "shell_bending_moment",
                                                   "shell_shear_force",
                                                   "shell_normal_force",
                                                   "thickness",
                                                   "shell_element_variables",
                                                   "strain_inner",
                                                   "strain_outer",
                                                   "internal_energy",
                                                   "beam_axial_force",
                                                   "beam_shear_force",
                                                   "beam_bending_moment",
                                                   "beam_torsion_moment",
                                                   "beam_axial_stress",
                                                   "beam_shear_stress",
                                                   "beam_axial_strain",
                                                   "sph_radius",
                                                   "sph_pressure",
                                                   "density",
                                                   "sph_neighbors",
                                                   "strain_rate",
                                                   "mass"};
        return v;
    }

    bool Has(int Family, std::initializer_list<const char*> Names) const {
        for (const char* n : Names)
            if (mSeen[Family].count(n))
                return true;
        return false;
    }

    std::size_t SeenPoints(int Family, const char* pName, std::size_t Default) const {
        const auto it = mSeen[Family].find(pName);
        return it == mSeen[Family].end() ? Default : it->second[0];
    }
    std::size_t SeenComps(int Family, const char* pName) const {
        const auto it = mSeen[Family].find(pName);
        return it == mSeen[Family].end() ? 0 : it->second[1];
    }

    void ElementLayouts() {
        const D3wModel& m = mModel;
        mStress = Has(kD3wSolid, {"stress"}) || Has(kD3wShell, {"stress"});
        mPstrain = Has(kD3wSolid, {"effective_plastic_strain"}) ||
                   Has(kD3wShell, {"effective_plastic_strain"});
        mPlasticTensor =
            Has(kD3wSolid, {"plastic_strain_tensor"}) || Has(kD3wShell, {"plastic_strain_tensor"});
        mThermalTensor =
            Has(kD3wSolid, {"thermal_strain_tensor"}) || Has(kD3wShell, {"thermal_strain_tensor"});
        mStrain = Has(kD3wSolid, {"strain"}) || Has(kD3wShell, {"strain_inner", "strain_outer"});
        mForces = false;
        for (const auto& f : d3w_shell_forces())
            mForces = mForces || Has(kD3wShell, {f.first});
        mExtra = Has(kD3wShell, {"thickness", "shell_element_variables", "internal_energy"});

        mSolid = D3wLayout{};
        if (m.mCounts[kD3wSolid]) {
            bool eight = false;
            for (const auto& kv : mSeen[kD3wSolid])
                eight = eight || kv.second[0] == 8;
            mSolid.mLayers = eight ? 8 : 1;
            const std::size_t hist = SeenComps(kD3wSolid, "history_variables");
            if (mStress)
                mSolid.mLayer.push_back({"stress", 6});
            if (mPstrain)
                mSolid.mLayer.push_back({"effective_plastic_strain", 1});
            if (mPlasticTensor)
                mSolid.mLayer.push_back({"plastic_strain_tensor", 6});
            if (mThermalTensor)
                mSolid.mLayer.push_back({"thermal_strain_tensor", 6});
            if (hist)
                mSolid.mLayer.push_back({"history_variables", hist});
            if (mStrain)
                mSolid.mLayer.push_back({"strain", 6});
        }
        mNeiph = 0;
        for (const auto& slot : mSolid.mLayer)
            if (slot.first != "stress" && slot.first != "effective_plastic_strain")
                mNeiph += slot.second;
        mNt3d = SeenComps(kD3wSolid, "thermal_variables");

        mShell = D3wLayout{};
        mNeips = 0;
        if (m.mCounts[kD3wShell]) {
            std::size_t layers = 0;
            for (const char* n : {"stress", "effective_plastic_strain", "history_variables"})
                layers = std::max(layers, SeenPoints(kD3wShell, n, 0));
            mShell.mLayers = layers ? layers : 1;
            const std::size_t hist = SeenComps(kD3wShell, "history_variables");
            if (mStress)
                mShell.mLayer.push_back({"stress", 6});
            if (mPstrain)
                mShell.mLayer.push_back({"effective_plastic_strain", 1});
            if (hist)
                mShell.mLayer.push_back({"history_variables", hist});
            mNeips = hist;
            if (mForces)
                for (const auto& f : d3w_shell_forces())
                    mShell.mTail.push_back({f.first, f.second, 1});
            if (mExtra)
                for (const auto& f : d3w_shell_extra())
                    mShell.mTail.push_back({f.first, f.second, 1});
            if (mStrain) {
                mShell.mTail.push_back({"strain_inner", 6, 1});
                mShell.mTail.push_back({"strain_outer", 6, 1});
            }
            if (mExtra)
                mShell.mTail.push_back({"internal_energy", 1, 1});
            if (mPlasticTensor)
                mShell.mTail.push_back(
                    {"plastic_strain_tensor", 6,
                     SeenPoints(kD3wShell, "plastic_strain_tensor", mShell.mLayers)});
            if (mThermalTensor)
                mShell.mTail.push_back({"thermal_strain_tensor", 6, 1});
        }

        mBeam = false;
        mBeamLayers = 0;
        mNeipb = 0;
        if (m.mCounts[kD3wBeam]) {
            std::size_t layers = 0;
            bool any_layer = false;
            for (const auto& l : d3w_beam_layer())
                if (mSeen[kD3wBeam].count(l.first)) {
                    any_layer = true;
                    layers = std::max(layers, SeenPoints(kD3wBeam, l.first, 0));
                }
            const std::size_t hist = SeenComps(kD3wBeam, "history_variables");
            bool core = false;
            for (const auto& c : d3w_beam_core())
                core = core || Has(kD3wBeam, {c.first});
            if (hist) {
                any_layer = true;
                const std::size_t hp = SeenPoints(kD3wBeam, "history_variables", 0);
                layers = std::max(layers, hp > 3 ? hp - 3 : std::size_t{0});
            }
            if (core || any_layer || hist) {
                mBeam = true;
                mBeamLayers = layers;
                mNeipb = hist;
            }
        }
        mSphFlags.assign(11, 0);
        mNsphVars = 0;
        if (m.mCounts[kD3wSph]) {
            std::vector<std::int64_t> flags(11, 0);
            flags[0] = 11;
            for (const D3wSphVar& var : d3w_sph_vars())
                if (mSeen[kD3wSph].count(var.mName))
                    flags[var.mFlag] = static_cast<std::int64_t>(var.mComps);
            if (mSeen[kD3wSph].count("strain_rate"))
                flags[8] = 12;
            else if (mSeen[kD3wSph].count("strain"))
                flags[8] = 6;
            flags[9] = mSeen[kD3wSph].count("mass") ? 1 : 0;
            flags[10] = static_cast<std::int64_t>(SeenComps(kD3wSph, "history_variables"));
            mSphFlags = flags;
            std::int64_t vars = 1;
            for (std::size_t i = 1; i < 8; ++i)
                vars += flags[i];
            vars += std::abs(flags[8]) + flags[9] + flags[10];
            mNsphVars = static_cast<std::size_t>(vars);
        }
        mNv3d = m.mCounts[kD3wSolid] ? mSolid.Width() : 0;
        mNv2d = m.mCounts[kD3wShell] ? mShell.Width() : 0;
        mNv1d = mBeam ? 6 + 5 * mBeamLayers + (3 + mBeamLayers) * mNeipb : 0;
    }

    void NodalFlags() {
        mDisp = mNode.count("displacement") != 0;
        mVel = mNode.count("velocity") != 0;
        mAcc = mNode.count("acceleration") != 0;
        mIt0 = 0;
        if (mNode.count("temperature") && mNode["temperature"] > 1)
            mIt0 = 3;
        else if (mNode.count("heat_flux"))
            mIt0 = 2;
        else if (mNode.count("temperature"))
            mIt0 = 1;
        mIt1 = mNode.count("mass_scaling") ? 1 : 0;
        mTgrad = mNode.count("temperature_gradient") != 0;
        mResid = mNode.count("residual_forces") || mNode.count("residual_moments");
    }

    std::size_t NodeWidth() const {
        std::size_t w = 3 * static_cast<std::size_t>(mDisp);
        if (mIt0) {
            w += mIt0 == 1 ? 1 : 4;
            if (mIt0 == 3)
                w += 2;
        }
        if (mIt1)
            w += 1;
        if (mTgrad)
            w += 1;
        if (mResid)
            w += 6;
        return w + 3 * static_cast<std::size_t>(mVel) + 3 * static_cast<std::size_t>(mAcc);
    }

    std::size_t StateWords() const {
        const D3wModel& m = mModel;
        std::size_t n = 1 + mNglbv + m.mN * NodeWidth();
        n += m.mCounts[kD3wSolid] * (mNv3d + mNt3d);
        n += m.mCounts[kD3wBeam] * mNv1d + m.mCounts[kD3wShell] * mNv2d;
        if (mAlive)
            n += m.mCounts[0] + m.mCounts[1] + m.mCounts[2];
        else if (mAliveNodes)
            n += m.mN;
        n += m.mCounts[kD3wSph] * mNsphVars + m.mNrigid * 24;
        return n;
    }

    std::int64_t MaxInt() const {
        const auto layers = static_cast<std::int64_t>(
            std::max<std::size_t>(mModel.mCounts[kD3wShell] ? mShell.mLayers : 1, 1));
        if (mAlive)
            return -(10000 + layers);
        if (mAliveNodes)
            return -layers;
        return layers;
    }

    std::int64_t Idtdt() const {
        return 100000 + 10000 * static_cast<std::int64_t>(mStrain) +
               1000 * static_cast<std::int64_t>(mThermalTensor) +
               100 * static_cast<std::int64_t>(mPlasticTensor) +
               10 * static_cast<std::int64_t>(mResid) + static_cast<std::int64_t>(mTgrad);
    }
    std::int64_t It() const { return 10 * mIt1 + mIt0; }
    std::size_t ExtraWords() const {
        return (mNt3d || mNeipb || !mModel.mExtra20.empty() || !mModel.mExtra27.empty()) ? 16 : 0;
    }

    const D3wModel& mModel;
    int mWs;
    std::map<std::string, std::size_t> mNode;
    std::array<D3wSeen, 4> mSeen;
    std::map<std::string, std::size_t> mRigid;  ///< reader name -> components, over all steps
    std::array<std::size_t, 3> mNummat{};
    std::size_t mPartCount = 0;
    std::vector<std::int64_t> mSphFlags = std::vector<std::int64_t>(11, 0);
    std::size_t mNsphVars = 0;
    bool mAlive = false, mAliveNodes = false, mGlobals = false, mPartGlobals = false;
    bool mStress = false, mPstrain = false, mPlasticTensor = false, mThermalTensor = false;
    bool mStrain = false, mForces = false, mExtra = false, mBeam = false;
    bool mDisp = false, mVel = false, mAcc = false, mTgrad = false, mResid = false;
    std::int64_t mIt0 = 0, mIt1 = 0;
    D3wLayout mSolid, mShell;
    std::size_t mNeiph = 0, mNeips = 0, mNeipb = 0, mNt3d = 0, mBeamLayers = 0;
    std::size_t mNglbv = 0, mNv3d = 0, mNv2d = 0, mNv1d = 0, mWords = 0;
};

std::string d3w_geometry(const D3wPlan& rPlan) {
    const D3wModel& m = rPlan.mModel;
    const int ws = rPlan.mWs;
    std::vector<std::int64_t> words(64 + rPlan.ExtraWords(), 0);
    words[11] = 1;
    words[15] = m.mNrigid ? 8 : 4;
    words[16] = static_cast<std::int64_t>(m.mN);
    words[17] = 6;
    words[18] = static_cast<std::int64_t>(rPlan.mNglbv);
    words[19] = rPlan.It();
    words[20] = rPlan.mDisp;
    words[21] = rPlan.mVel;
    words[22] = rPlan.mAcc;
    words[23] = m.mExtra.empty() ? static_cast<std::int64_t>(m.mCounts[kD3wSolid])
                                 : -static_cast<std::int64_t>(m.mCounts[kD3wSolid]);
    words[24] = static_cast<std::int64_t>(rPlan.mNummat[kD3wSolid]);
    words[27] = static_cast<std::int64_t>(rPlan.mNv3d);
    words[28] = static_cast<std::int64_t>(m.mCounts[kD3wBeam]);
    words[29] = static_cast<std::int64_t>(m.mNummat[kD3wBeam]);
    words[30] = static_cast<std::int64_t>(rPlan.mNv1d);
    words[31] = static_cast<std::int64_t>(m.mCounts[kD3wShell]);
    words[32] = static_cast<std::int64_t>(m.mNummat[kD3wShell]);
    words[33] = static_cast<std::int64_t>(rPlan.mNv2d);
    words[34] = static_cast<std::int64_t>(rPlan.mNeiph);
    words[35] = static_cast<std::int64_t>(rPlan.mNeips);
    words[36] = rPlan.MaxInt();
    words[37] = static_cast<std::int64_t>(m.mCounts[kD3wSph]);
    words[38] = m.mCounts[kD3wSph] ? 1 : 0;
    words[43] = rPlan.mStress ? 1000 : 0;
    words[44] = rPlan.mPstrain ? 1000 : 0;
    words[45] = rPlan.mForces ? 1000 : 0;
    words[46] = rPlan.mExtra ? 1000 : 0;
    words[51] = static_cast<std::int64_t>(m.mPartIds.size());
    words[55] = static_cast<std::int64_t>(m.mExtra8.size());
    words[56] = rPlan.Idtdt();
    if (rPlan.ExtraWords()) {
        words[57] = static_cast<std::int64_t>(rPlan.ExtraWords());
        words[64] = static_cast<std::int64_t>(m.mExtra20.size());
        words[65] = static_cast<std::int64_t>(rPlan.mNt3d);
        words[66] = static_cast<std::int64_t>(m.mExtra27.size());
        words[67] = static_cast<std::int64_t>(rPlan.mNeipb);
    }
    const std::size_t ids = m.mCounts[0] + m.mCounts[1] + m.mCounts[2];
    words[39] = static_cast<std::int64_t>(16 + m.mN + ids + 3 * m.mPartIds.size());
    std::string out;
    for (const std::int64_t w : words)
        d3w_i32(out, w, ws);
    const std::string title = "meshio++" + std::string(10 * static_cast<std::size_t>(ws) - 8, ' ');
    out.replace(0, title.size(), title);
    std::string version;
    d3w_f32(version, kD3wVersion, ws);
    out.replace(14 * static_cast<std::size_t>(ws), static_cast<std::size_t>(ws), version);
    if (m.mCounts[kD3wSph])
        for (const std::int64_t f : rPlan.mSphFlags)
            d3w_i32(out, f, ws);
    for (const double v : m.mPoints)
        d3w_f32(out, v, ws);
    for (std::size_t e = 0; e < m.mCounts[kD3wSolid]; ++e) {
        for (const std::int64_t v : m.mRows[kD3wSolid][e])
            d3w_i32(out, v + 1, ws);
        d3w_i32(out, m.mPartIndex.at(m.mPartsOf[kD3wSolid][e]), ws);
    }
    if (!m.mExtra.empty()) {
        std::vector<std::array<std::int64_t, 2>> extra(m.mCounts[kD3wSolid], {0, 0});
        for (const auto& [e, nodes] : m.mExtra)
            extra[e] = {nodes[0] + 1, nodes[1] + 1};
        for (const auto& row : extra) {
            d3w_i32(out, row[0], ws);
            d3w_i32(out, row[1], ws);
        }
    }
    for (std::size_t e = 0; e < m.mCounts[kD3wBeam]; ++e) {
        d3w_i32(out, m.mRows[kD3wBeam][e][0] + 1, ws);
        d3w_i32(out, m.mRows[kD3wBeam][e][1] + 1, ws);
        for (int k = 0; k < 3; ++k)
            d3w_i32(out, 0, ws);
        d3w_i32(out, m.mPartIndex.at(m.mPartsOf[kD3wBeam][e]), ws);
    }
    for (std::size_t e = 0; e < m.mCounts[kD3wShell]; ++e) {
        for (std::size_t k = 0; k < 4; ++k)
            d3w_i32(out, m.mRows[kD3wShell][e][k] + 1, ws);
        d3w_i32(out, m.mPartIndex.at(m.mPartsOf[kD3wShell][e]), ws);
    }
    const auto nparts = static_cast<std::int64_t>(m.mPartIds.size());
    const auto solids_at = static_cast<std::int64_t>(17 + m.mN);
    const auto beams_at = solids_at + static_cast<std::int64_t>(m.mCounts[kD3wSolid]);
    const auto shells_at = beams_at + static_cast<std::int64_t>(m.mCounts[kD3wBeam]);
    const auto parts_at = shells_at + static_cast<std::int64_t>(m.mCounts[kD3wShell]);
    for (const std::int64_t v :
         {std::int64_t{-17}, solids_at, beams_at, shells_at, parts_at,
          static_cast<std::int64_t>(m.mN), static_cast<std::int64_t>(m.mCounts[kD3wSolid]),
          static_cast<std::int64_t>(m.mCounts[kD3wBeam]),
          static_cast<std::int64_t>(m.mCounts[kD3wShell]), std::int64_t{0}, parts_at,
          parts_at + nparts, parts_at + 2 * nparts, nparts, static_cast<std::int64_t>(m.mNrigid),
          nparts})
        d3w_i32(out, v, ws);
    for (const std::int64_t v : m.mNid)
        d3w_i32(out, v, ws);
    for (const int family : {kD3wSolid, kD3wBeam, kD3wShell})
        for (const std::int64_t v : m.mEid[family])
            d3w_i32(out, v, ws);
    for (int copy = 0; copy < 2; ++copy)
        for (const std::int64_t v : m.mPartIds)
            d3w_i32(out, v, ws);
    for (std::int64_t k = 1; k <= nparts; ++k)
        d3w_i32(out, k, ws);
    // rigid bodies (part, no node lists), then the SPH particles' (node, part)
    if (m.mNrigid) {
        d3w_i32(out, static_cast<std::int64_t>(m.mNrigid), ws);
        for (const std::int64_t part : m.mRigidParts) {
            d3w_i32(out, m.mPartIndex.at(part), ws);
            d3w_i32(out, 0, ws);
            d3w_i32(out, 0, ws);
        }
    }
    for (std::size_t e = 0; e < m.mCounts[kD3wSph]; ++e) {
        d3w_i32(out, m.mRows[kD3wSph][e][0] + 1, ws);
        d3w_i32(out, m.mPartIndex.at(m.mPartsOf[kD3wSph][e]), ws);
    }
    // the higher-order tables: eight-node shells, then 20- and 27-node solids
    for (const auto* table : {&m.mExtra8, &m.mExtra20, &m.mExtra27})
        for (const auto& [e, nodes] : *table) {
            d3w_i32(out, static_cast<std::int64_t>(e) + 1, ws);
            for (const std::int64_t v : nodes)
                d3w_i32(out, v + 1, ws);
        }
    d3w_f32(out, kD3EofMarker, ws);
    d3w_i32(out, 90001, ws);
    d3w_i32(out, static_cast<std::int64_t>(m.mPartIds.size()), ws);
    for (const std::int64_t pid : m.mPartIds) {
        const auto it = m.mTitles.find(pid);
        std::string text = it != m.mTitles.end() ? it->second : "Part " + std::to_string(pid);
        text = text.substr(0, 72);
        d3w_i32(out, pid, ws);
        out += text;
        out.append(72 - text.size(), ' ');
    }
    d3w_f32(out, kD3EofMarker, ws);
    return out;
}

/// A point array as (n, comps) doubles: missing components and arrays are NaN.
std::vector<double> d3w_node_array(const Mesh& rMesh, const char* pName, std::size_t N,
                                   std::size_t Comps = 3) {
    std::vector<double> out(N * Comps, kD3Nan);
    if (!rMesh.HasPointData(pName))
        return out;
    const NDArray& a = rMesh.PointData(pName);
    const std::size_t width = a.Size() / N;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t k = 0; k < Comps && k < width; ++k)
            out[i * Comps + k] = detail::read_double(a, i * width + k);
    return out;
}

/// A family's element variables as (cells, width), NaN where absent.
std::vector<double> d3w_family_matrix(const D3wPlan& rPlan, const Mesh& rMesh, int Family,
                                      const D3wLayout& rLayout) {
    const D3wModel& m = rPlan.mModel;
    const std::size_t width = rLayout.Width();
    std::vector<double> out(m.mCounts[Family] * width, kD3Nan);
    const std::size_t layer_width = rLayout.LayerWidth();
    std::size_t offset = 0;
    for (const auto& [name, comps] : rLayout.mLayer) {
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, name, present);
        if (present)
            for (const D3wRange& range : m.mRanges[Family]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, comps);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t layer = 0; layer < std::min(a.mPoints, rLayout.mLayers);
                         ++layer)
                        for (std::size_t c = 0; c < k; ++c)
                            out[(range.mFirst + e) * width + layer * layer_width + offset + c] =
                                a.mV[(e * a.mPoints + layer) * a.mComps + c];
            }
        offset += comps;
    }
    std::size_t base = rLayout.mLayers * layer_width;
    for (const auto& tail : rLayout.mTail) {
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, tail.mName, present);
        if (present)
            for (const D3wRange& range : m.mRanges[Family]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, tail.mComps);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t point = 0; point < std::min(a.mPoints, tail.mCount); ++point)
                        for (std::size_t c = 0; c < k; ++c)
                            out[(range.mFirst + e) * width + base + point * tail.mComps + c] =
                                a.mV[(e * a.mPoints + point) * a.mComps + c];
            }
        base += tail.mComps * tail.mCount;
    }
    return out;
}

/// The beams' columns: six resultants, then per layer the axial stress, shear stress (2),
/// plastic strain and axial strain, then the history values.
std::vector<double> d3w_beam_matrix(const D3wPlan& rPlan, const Mesh& rMesh) {
    const D3wModel& m = rPlan.mModel;
    const std::size_t nl = rPlan.mBeamLayers, nh = rPlan.mNeipb, width = rPlan.mNv1d;
    std::vector<double> out(m.mCounts[kD3wBeam] * width, kD3Nan);
    std::size_t at = 0;
    for (const auto& [name, comps] : d3w_beam_core()) {
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, name, present);
        if (present)
            for (const D3wRange& range : m.mRanges[kD3wBeam]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, comps);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t c = 0; c < k; ++c)
                        out[(range.mFirst + e) * width + at + c] =
                            a.mV[e * a.mPoints * a.mComps + c];
            }
        at += comps;
    }
    std::size_t offset = 0;
    for (const auto& [name, comps] : d3w_beam_layer()) {
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, name, present);
        if (present)
            for (const D3wRange& range : m.mRanges[kD3wBeam]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, comps);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t layer = 0; layer < std::min(a.mPoints, nl); ++layer)
                        for (std::size_t c = 0; c < k; ++c)
                            out[(range.mFirst + e) * width + 6 + layer * 5 + offset + c] =
                                a.mV[(e * a.mPoints + layer) * a.mComps + c];
            }
        offset += comps;
    }
    if (nh) {
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, "history_variables", present);
        if (present)
            for (const D3wRange& range : m.mRanges[kD3wBeam]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, nh);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t point = 0; point < std::min(a.mPoints, 3 + nl); ++point)
                        for (std::size_t c = 0; c < k; ++c)
                            out[(range.mFirst + e) * width + 6 + 5 * nl + point * nh + c] =
                                a.mV[(e * a.mPoints + point) * a.mComps + c];
            }
    }
    return out;
}

/// The SPH particles' rows: the material (negative once deleted), then the variables the flag
/// block names.
std::vector<double> d3w_sph_matrix(const D3wPlan& rPlan, const Mesh& rMesh) {
    const D3wModel& m = rPlan.mModel;
    const std::size_t width = rPlan.mNsphVars;
    std::vector<double> out(m.mCounts[kD3wSph] * width, kD3Nan);
    for (std::size_t e = 0; e < m.mCounts[kD3wSph]; ++e)
        out[e * width] = static_cast<double>(m.mPartIndex.at(m.mPartsOf[kD3wSph][e]));
    if (rMesh.HasCellData("lsdyna:alive"))
        for (const D3wRange& range : m.mRanges[kD3wSph]) {
            if (range.mBlock >= rMesh.CellDataNumBlocks("lsdyna:alive"))
                continue;
            const NDArray& a = rMesh.CellData("lsdyna:alive", range.mBlock);
            for (std::size_t e = 0; e < range.mCount && e < a.Size(); ++e)
                if (detail::read_double(a, e) == 0.0)
                    out[(range.mFirst + e) * width] = -out[(range.mFirst + e) * width];
        }
    const auto& flags = rPlan.mSphFlags;
    const auto strain = static_cast<std::size_t>(std::abs(flags[8]));
    const std::vector<std::pair<const char*, std::size_t>> columns = {
        {"sph_radius", static_cast<std::size_t>(flags[1])},
        {"sph_pressure", static_cast<std::size_t>(flags[2])},
        {"stress", static_cast<std::size_t>(flags[3])},
        {"effective_plastic_strain", static_cast<std::size_t>(flags[4])},
        {"density", static_cast<std::size_t>(flags[5])},
        {"internal_energy", static_cast<std::size_t>(flags[6])},
        {"sph_neighbors", static_cast<std::size_t>(flags[7])},
        {"strain", std::min<std::size_t>(strain, 6)},
        {"strain_rate", strain - std::min<std::size_t>(strain, 6)},
        {"mass", static_cast<std::size_t>(flags[9])},
        {"history_variables", static_cast<std::size_t>(flags[10])}};
    std::size_t at = 1;
    for (const auto& [name, w] : columns) {
        if (w > 0) {
            bool present = false;
            const auto arrays = d3w_cell_array(rMesh, name, present);
            if (present)
                for (const D3wRange& range : m.mRanges[kD3wSph]) {
                    const D3wBlockArray& a = arrays[range.mBlock];
                    if (a.mPoints == 0 || range.mCount == 0)
                        continue;
                    const std::size_t k = std::min(a.mComps, w);
                    for (std::size_t e = 0; e < range.mCount; ++e)
                        for (std::size_t c = 0; c < k; ++c)
                            out[(range.mFirst + e) * width + at + c] =
                                a.mV[e * a.mPoints * a.mComps + c];
                }
        }
        at += w;
    }
    return out;
}

std::string d3w_state(const D3wPlan& rPlan, const Mesh& rMesh, double Time) {
    const D3wModel& m = rPlan.mModel;
    const std::size_t n = m.mN;
    const int ws = rPlan.mWs;
    std::string out;
    d3w_f32(out, Time, ws);
    if (rPlan.mNglbv) {
        std::vector<double> g(rPlan.mNglbv, 0.0);
        std::size_t at = 0;
        for (const auto& [name, width] : d3w_globals()) {
            if (rMesh.HasFieldData(name)) {
                const NDArray& v = rMesh.FieldData(name);
                for (std::size_t k = 0; k < width && k < v.Size(); ++k)
                    g[at + k] = detail::read_double(v, k);
            }
            at += width;
        }
        if (rPlan.mPartGlobals) {
            const std::size_t nparts = rPlan.mPartCount;
            for (const auto& [name, w] : d3w_part_globals()) {
                if (rMesh.HasFieldData(name) && rMesh.FieldData(name).Size() == w * nparts)
                    for (std::size_t k = 0; k < w * nparts; ++k)
                        g[at + k] = detail::read_double(rMesh.FieldData(name), k);
                at += w * nparts;
            }
        }
        for (const double v : g)
            d3w_f32(out, v, ws);
    }
    // nodal variables, in the reader's order
    if (rPlan.mDisp) {
        const std::vector<double> d = d3w_node_array(rMesh, "displacement", n);
        for (std::size_t i = 0; i < 3 * n; ++i)
            d3w_f32(out, m.mPoints[i] + (std::isnan(d[i]) ? 0.0 : d[i]), ws);
    }
    if (rPlan.mIt0) {
        for (const double v : d3w_node_array(rMesh, "temperature", n, rPlan.mIt0 == 3 ? 3 : 1))
            d3w_f32(out, v, ws);
        if (rPlan.mIt0 >= 2)
            for (const double v : d3w_node_array(rMesh, "heat_flux", n))
                d3w_f32(out, v, ws);
    }
    if (rPlan.mIt1)
        for (const double v : d3w_node_array(rMesh, "mass_scaling", n, 1))
            d3w_f32(out, v, ws);
    if (rPlan.mTgrad)
        for (const double v : d3w_node_array(rMesh, "temperature_gradient", n, 1))
            d3w_f32(out, v, ws);
    if (rPlan.mResid) {
        for (const double v : d3w_node_array(rMesh, "residual_forces", n))
            d3w_f32(out, v, ws);
        for (const double v : d3w_node_array(rMesh, "residual_moments", n))
            d3w_f32(out, v, ws);
    }
    if (rPlan.mVel)
        for (const double v : d3w_node_array(rMesh, "velocity", n))
            d3w_f32(out, v, ws);
    if (rPlan.mAcc)
        for (const double v : d3w_node_array(rMesh, "acceleration", n))
            d3w_f32(out, v, ws);
    if (rPlan.mNt3d) {
        std::vector<double> thermal(m.mCounts[kD3wSolid] * rPlan.mNt3d, kD3Nan);
        bool present = false;
        const auto arrays = d3w_cell_array(rMesh, "thermal_variables", present);
        if (present)
            for (const D3wRange& range : m.mRanges[kD3wSolid]) {
                const D3wBlockArray& a = arrays[range.mBlock];
                if (a.mPoints == 0 || range.mCount == 0)
                    continue;
                const std::size_t k = std::min(a.mComps, rPlan.mNt3d);
                for (std::size_t e = 0; e < range.mCount; ++e)
                    for (std::size_t c = 0; c < k; ++c)
                        thermal[(range.mFirst + e) * rPlan.mNt3d + c] =
                            a.mV[e * a.mPoints * a.mComps + c];
            }
        for (const double v : thermal)
            d3w_f32(out, v, ws);
    }
    if (rPlan.mNv3d)
        for (const double v : d3w_family_matrix(rPlan, rMesh, kD3wSolid, rPlan.mSolid))
            d3w_f32(out, v, ws);
    if (rPlan.mNv1d)
        for (const double v : d3w_beam_matrix(rPlan, rMesh))
            d3w_f32(out, v, ws);
    if (rPlan.mNv2d)
        for (const double v : d3w_family_matrix(rPlan, rMesh, kD3wShell, rPlan.mShell))
            d3w_f32(out, v, ws);
    if (rPlan.mAlive) {
        const bool has = rMesh.HasCellData("lsdyna:alive");
        for (const int family : {kD3wSolid, kD3wShell, kD3wBeam}) {
            std::vector<double> flags(m.mCounts[family], 1.0);
            if (has)
                for (const D3wRange& range : m.mRanges[family]) {
                    if (range.mBlock >= rMesh.CellDataNumBlocks("lsdyna:alive"))
                        continue;
                    const NDArray& a = rMesh.CellData("lsdyna:alive", range.mBlock);
                    for (std::size_t e = 0; e < range.mCount && e < a.Size(); ++e)
                        flags[range.mFirst + e] = detail::read_double(a, e) != 0.0 ? 1.0 : 0.0;
                }
            for (const double v : flags)
                d3w_f32(out, v, ws);
        }
    } else if (rPlan.mAliveNodes) {
        std::vector<double> flags(n, 1.0);
        if (rMesh.HasPointData("lsdyna:alive")) {
            const NDArray& a = rMesh.PointData("lsdyna:alive");
            for (std::size_t i = 0; i < n && i < a.Size(); ++i)
                flags[i] = detail::read_double(a, i) != 0.0 ? 1.0 : 0.0;
        }
        for (const double v : flags)
            d3w_f32(out, v, ws);
    }
    if (m.mCounts[kD3wSph])
        for (const double v : d3w_sph_matrix(rPlan, rMesh))
            d3w_f32(out, v, ws);
    if (m.mNrigid) {
        std::vector<double> rigid(m.mNrigid * 24, kD3Nan);
        std::size_t at = 0;
        for (const auto& [name, width] : d3w_rigid_fields()) {
            const std::string key = std::string("lsdyna:rigid_body_") + name;
            if (rMesh.HasFieldData(key) && rMesh.FieldData(key).Size() == m.mNrigid * width)
                for (std::size_t r = 0; r < m.mNrigid; ++r)
                    for (std::size_t c = 0; c < width; ++c)
                        rigid[r * 24 + at + c] =
                            detail::read_double(rMesh.FieldData(key), r * width + c);
            at += width;
        }
        for (const double v : rigid)
            d3w_f32(out, v, ws);
    }
    return out;
}

/// The files of the family, each with the geometry in the first and the end mark last.
std::vector<std::string> d3w_assemble(const D3wPlan& rPlan,
                                      const std::vector<std::pair<double, const Mesh*>>& rSteps,
                                      std::size_t MaxBytes) {
    std::string mark;
    d3w_f32(mark, kD3EofMarker, rPlan.mWs);
    std::vector<std::string> files;
    std::string current = d3w_geometry(rPlan);
    bool any = false;
    for (const auto& [time, pMesh] : rSteps) {
        const std::string state = d3w_state(rPlan, *pMesh, time);
        if (any && current.size() + state.size() + mark.size() > MaxBytes) {
            files.push_back(std::move(current) + mark);
            current.clear();
            any = false;
        }
        current += state;
        any = true;
    }
    files.push_back(std::move(current) + mark);
    return files;
}

void d3w_store(const std::string& rPath, const std::vector<std::string>& rFiles) {
    const std::vector<std::string> stale = d3_family(rPath);
    for (std::size_t k = 1; k < stale.size(); ++k) {
        std::error_code ec;
        fs::remove(stale[k], ec);
    }
    for (std::size_t k = 0; k < rFiles.size(); ++k) {
        std::string suffix = std::to_string(k);
        if (suffix.size() < 2)
            suffix.insert(0, "0");
        const std::string target = k == 0 ? rPath : rPath + suffix;
        auto out = detail::make_classic_ofstream(target, std::ios::binary);
        if (!out)
            throw WriteError("Could not open file for writing: " + target);
        out.write(rFiles[k].data(), static_cast<std::streamsize>(rFiles[k].size()));
        out.close();
        if (!out)
            d3w_fail("failed writing " + target);
    }
}

}  // namespace

void write_lsdyna_d3plot(const std::string& rPath, const Mesh& rMesh, std::size_t MaxFileBytes,
                         bool DoublePrecision) {
    const D3wModel model(rMesh);
    const D3wPlan plan(model, {&rMesh}, DoublePrecision);
    d3w_store(rPath, d3w_assemble(plan, {{d3w_time_of(rMesh), &rMesh}}, MaxFileBytes));
}

struct LsdynaD3plotSeriesWriter::Impl {
    std::string mPath;
    std::size_t mMaxBytes = kD3wDefaultFileBytes;
    bool mDouble = false;
    std::unique_ptr<D3wModel> mModel;
    std::vector<double> mTimes;
    std::vector<std::unique_ptr<Mesh>> mMeshes;
    bool mDone = false;
};

LsdynaD3plotSeriesWriter::LsdynaD3plotSeriesWriter(const std::string& rPath,
                                                   std::size_t MaxFileBytes, bool DoublePrecision)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
    mpImpl->mMaxBytes = MaxFileBytes;
    mpImpl->mDouble = DoublePrecision;
}

LsdynaD3plotSeriesWriter::~LsdynaD3plotSeriesWriter() = default;
LsdynaD3plotSeriesWriter::LsdynaD3plotSeriesWriter(LsdynaD3plotSeriesWriter&&) noexcept = default;
LsdynaD3plotSeriesWriter& LsdynaD3plotSeriesWriter::operator=(LsdynaD3plotSeriesWriter&&) noexcept =
    default;

void LsdynaD3plotSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        d3w_fail("series: the writer was moved from");
    Impl& s = *mpImpl;
    auto model = std::make_unique<D3wModel>(rMesh);
    if (!s.mModel)
        s.mModel = std::move(model);
    else if (model->mSignature != s.mModel->mSignature)
        d3w_fail("series: step " + std::to_string(s.mMeshes.size()) +
                 "'s points or cells differ from the first step's; a family holds one mesh -- "
                 "write one family per step with '{step}'");
    s.mTimes.push_back(Time);
    s.mMeshes.push_back(std::make_unique<Mesh>(detail::clone_mesh(rMesh)));
}

std::size_t LsdynaD3plotSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mMeshes.size() : 0;
}

void LsdynaD3plotSeriesWriter::Finalize() {
    if (!mpImpl || mpImpl->mDone)
        return;
    Impl& s = *mpImpl;
    if (s.mMeshes.empty())
        d3w_fail("series: no step was written");
    s.mDone = true;
    std::vector<const Mesh*> meshes;
    std::vector<std::pair<double, const Mesh*>> steps;
    for (std::size_t i = 0; i < s.mMeshes.size(); ++i) {
        meshes.push_back(s.mMeshes[i].get());
        steps.push_back({s.mTimes[i], s.mMeshes[i].get()});
    }
    const D3wPlan plan(*s.mModel, meshes, s.mDouble);
    d3w_store(s.mPath, d3w_assemble(plan, steps, s.mMaxBytes));
}

}  // namespace meshioplusplus
