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
#include <cstring>
#include <filesystem>
#include <ios>
#include <iterator>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/radioss_anim.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"

namespace meshioplusplus {

namespace {

constexpr std::int32_t kAnimMagic = 0x542C;

// Big-endian reads with a bounds check.
class AnimCursor {
public:
    AnimCursor(std::string_view rData, const std::string& rPath) : mData(rData), mPath(rPath) {}

    const char* Take(std::size_t N) {
        if (N > mData.size() - mPos)  // not mPos + N: that wraps for a huge N
            throw ReadError("Radioss animation: '" + mPath + "' is truncated (needs " +
                            std::to_string(N) + " more bytes at offset " + std::to_string(mPos) +
                            " of " + std::to_string(mData.size()) + ")");
        const char* p = mData.data() + mPos;
        mPos += N;
        return p;
    }

    std::vector<std::int64_t> Ints(std::int64_t N) {
        if (N < 0)
            throw ReadError("Radioss animation: '" + mPath + "' has a negative count");
        if (static_cast<std::uint64_t>(N) > (mData.size() - mPos) / 4)
            Take(mData.size() - mPos + 1);  // reports the truncation
        const unsigned char* p =
            reinterpret_cast<const unsigned char*>(Take(4 * static_cast<std::size_t>(N)));
        std::vector<std::int64_t> out(static_cast<std::size_t>(N));
        for (std::size_t k = 0; k < out.size(); ++k) {
            const std::uint32_t v =
                (std::uint32_t(p[4 * k]) << 24) | (std::uint32_t(p[4 * k + 1]) << 16) |
                (std::uint32_t(p[4 * k + 2]) << 8) | std::uint32_t(p[4 * k + 3]);
            out[k] = static_cast<std::int32_t>(v);
        }
        return out;
    }

    std::int64_t Int() { return Ints(1)[0]; }

    std::vector<double> Floats(std::int64_t N) {
        if (N < 0)
            throw ReadError("Radioss animation: '" + mPath + "' has a negative count");
        if (static_cast<std::uint64_t>(N) > (mData.size() - mPos) / 4)
            Take(mData.size() - mPos + 1);  // reports the truncation
        const unsigned char* p =
            reinterpret_cast<const unsigned char*>(Take(4 * static_cast<std::size_t>(N)));
        std::vector<double> out(static_cast<std::size_t>(N));
        for (std::size_t k = 0; k < out.size(); ++k) {
            const std::uint32_t v =
                (std::uint32_t(p[4 * k]) << 24) | (std::uint32_t(p[4 * k + 1]) << 16) |
                (std::uint32_t(p[4 * k + 2]) << 8) | std::uint32_t(p[4 * k + 3]);
            float f;
            std::memcpy(&f, &v, 4);
            out[k] = static_cast<double>(f);
        }
        return out;
    }

    std::vector<bool> Flags(std::int64_t N) {
        const char* p = Take(static_cast<std::size_t>(N < 0 ? 0 : N));
        std::vector<bool> out(static_cast<std::size_t>(N < 0 ? 0 : N));
        for (std::size_t k = 0; k < out.size(); ++k)
            out[k] = p[k] != 0;
        return out;
    }

    std::string Text(std::size_t N) {
        std::string s(Take(N), N);
        const std::size_t nul = s.find('\0');
        if (nul != std::string::npos)
            s.resize(nul);
        const std::size_t b = s.find_first_not_of(" \t\r\n");
        if (b == std::string::npos)
            return {};
        const std::size_t e = s.find_last_not_of(" \t\r\n");
        return s.substr(b, e - b + 1);
    }

    std::vector<std::string> Texts(std::int64_t Count, std::size_t N) {
        std::vector<std::string> out;
        for (std::int64_t k = 0; k < Count; ++k)
            out.push_back(Text(N));
        return out;
    }

private:
    std::string_view mData;
    const std::string& mPath;
    std::size_t mPos = 0;
};

// One element family (1-D, 2-D, 3-D or SPH) of the file.
struct AnimFamily {
    std::string mKey;
    std::size_t mNodesPer = 0;
    std::vector<std::int64_t> mConn;  // mNodesPer per element
    std::vector<bool> mAlive;
    std::vector<std::int64_t> mPartEnds;
    std::vector<std::string> mPartTexts;
    std::vector<std::pair<std::string, std::vector<double>>> mScalars;
    // name, width (3, 6 or 9), values (width per element)
    std::vector<std::tuple<std::string, std::size_t, std::vector<double>>> mTensors;
    std::vector<double> mMass;
    std::vector<std::int64_t> mIds;
    std::vector<std::vector<std::int64_t>> mHierarchy;  // subset, material, property per part
    bool mHasMass = false, mHasIds = false, mHasHierarchy = false;

    std::size_t Size() const { return mNodesPer ? mConn.size() / mNodesPer : 0; }
};

// The shared tail of the 3-D, 1-D and SPH sections, after their counts and
// connectivity: alive flags, parts, scalars, tensors, then (1-D) the skews,
// masses, ids and the part hierarchy.
void anim_section(AnimCursor& rC, const std::vector<std::int64_t>& rFlags, AnimFamily& rF,
                  std::int64_t NumParts, std::int64_t NumScalars, std::size_t TensorWidth,
                  std::int64_t NumTensors, bool Skews) {
    const std::int64_t n = static_cast<std::int64_t>(rF.Size());
    rF.mAlive = rC.Flags(n);
    rF.mPartEnds = rC.Ints(NumParts);
    rF.mPartTexts = rC.Texts(NumParts, 50);
    if (NumScalars) {
        const std::vector<std::string> names = rC.Texts(NumScalars, 81);
        const std::vector<double> values = rC.Floats(NumScalars * n);
        for (std::size_t k = 0; k < names.size(); ++k)
            rF.mScalars.emplace_back(
                names[k],
                std::vector<double>(values.begin() + static_cast<std::ptrdiff_t>(k) * n,
                                    values.begin() + static_cast<std::ptrdiff_t>(k + 1) * n));
    }
    if (NumTensors) {
        const std::vector<std::string> names = rC.Texts(NumTensors, 81);
        const std::int64_t w = static_cast<std::int64_t>(TensorWidth);
        const std::vector<double> values = rC.Floats(n * w * NumTensors);
        for (std::size_t k = 0; k < names.size(); ++k)
            rF.mTensors.emplace_back(
                names[k], TensorWidth,
                std::vector<double>(values.begin() + static_cast<std::ptrdiff_t>(k) * n * w,
                                    values.begin() + static_cast<std::ptrdiff_t>(k + 1) * n * w));
    }
    if (Skews)
        rC.Ints(n);
    if (rFlags[0] == 1) {
        rF.mMass = rC.Floats(n);
        rF.mHasMass = true;
    }
    if (rFlags[1] == 1) {
        rF.mIds = rC.Ints(n);
        rF.mHasIds = true;
    }
    if (rFlags[4]) {
        for (int k = 0; k < 3; ++k)
            rF.mHierarchy.push_back(rC.Ints(NumParts));
        rF.mHasHierarchy = true;
    }
}

// Per element, part k runs up to (not including) `rEnds[k]`.
template <class Fn>
void anim_for_parts(const AnimFamily& rF, Fn&& rFn) {
    std::int64_t start = 0;
    for (std::size_t k = 0; k < rF.mPartEnds.size(); ++k) {
        const std::int64_t end = std::max(start, rF.mPartEnds[k]);
        rFn(k, start, std::min<std::int64_t>(end, static_cast<std::int64_t>(rF.Size())));
        start = end;
    }
}

}  // namespace

bool is_radioss_anim_filename(const std::string& rPath) {
    const std::string name = std::filesystem::path(rPath).filename().string();
    std::size_t k = name.size();
    while (k > 0 && name[k - 1] >= '0' && name[k - 1] <= '9')
        --k;
    return name.size() - k >= 3 && k >= 2 && name[k - 1] == 'A';
}

Mesh read_radioss_anim(const std::string& rPath) {
    const detail::FileSource data_source =
        detail::open_source(rPath, "Radioss animation: cannot open " + rPath);
    const std::string_view data = data_source.View();
    AnimCursor c(data, rPath);
    const std::int64_t magic = c.Int();
    if (magic != kAnimMagic) {
        std::string hex;
        for (std::uint32_t v = static_cast<std::uint32_t>(magic); v || hex.empty(); v >>= 4)
            hex.insert(hex.begin(), "0123456789abcdef"[v & 0xf]);
        throw ReadError("Radioss animation: '" + rPath + "' has magic 0x" + hex +
                        "; only the current layout (0x542c) is read");
    }
    const double time = c.Floats(1)[0];
    c.Texts(3, 81);  // "Time=", "ModAnim", "Radioss Run="
    const std::vector<std::int64_t> flags = c.Ints(10);
    const std::vector<std::int64_t> counts = c.Ints(8);
    const std::int64_t nn = counts[0], nf = counts[1], np2 = counts[2], nfun = counts[3],
                       nefun = counts[4], nvec = counts[5], nten = counts[6], nskew = counts[7];
    // Every count sizes records of at least one byte per entry, so none can
    // exceed the file; this also keeps the products below inside int64.
    for (const std::int64_t n : counts)
        detail::checked_count(n, data.size(), "Radioss animation", "header");
    c.Take(static_cast<std::size_t>(2 * 6 * std::max<std::int64_t>(nskew, 0)));
    const std::vector<double> coords = c.Floats(3 * nn);

    std::vector<AnimFamily> families;
    AnimFamily two;
    two.mKey = "2d";
    two.mNodesPer = 4;
    two.mConn = c.Ints(4 * nf);
    two.mAlive = c.Flags(nf);
    if (np2) {
        two.mPartEnds = c.Ints(np2);
        two.mPartTexts = c.Texts(np2, 50);
    }
    c.Take(static_cast<std::size_t>(2 * 3 * nn));  // nodal normals
    std::vector<std::pair<std::string, std::vector<double>>> point_scalars;
    if (nfun + nefun) {
        const std::vector<std::string> names = c.Texts(nfun + nefun, 81);
        const std::vector<double> values = c.Floats(nn * nfun);
        for (std::int64_t k = 0; k < nfun; ++k)
            point_scalars.emplace_back(
                names[static_cast<std::size_t>(k)],
                std::vector<double>(values.begin() + k * nn, values.begin() + (k + 1) * nn));
        const std::vector<double> evalues = c.Floats(nf * nefun);
        for (std::int64_t k = 0; k < nefun; ++k)
            two.mScalars.emplace_back(
                names[static_cast<std::size_t>(nfun + k)],
                std::vector<double>(evalues.begin() + k * nf, evalues.begin() + (k + 1) * nf));
    }
    const std::vector<std::string> vnames = c.Texts(nvec, 81);
    const std::vector<double> vvalues = c.Floats(3 * nn * nvec);
    if (nten) {
        const std::vector<std::string> names = c.Texts(nten, 81);
        const std::vector<double> values = c.Floats(nf * 3 * nten);
        for (std::int64_t k = 0; k < nten; ++k)
            two.mTensors.emplace_back(names[static_cast<std::size_t>(k)], 3,
                                      std::vector<double>(values.begin() + k * nf * 3,
                                                          values.begin() + (k + 1) * nf * 3));
    }
    std::vector<double> point_mass;
    if (flags[0] == 1) {
        two.mMass = c.Floats(nf);
        two.mHasMass = true;
        point_mass = c.Floats(nn);
    }
    std::vector<std::int64_t> node_ids;
    if (flags[1]) {
        node_ids = c.Ints(nn);
        two.mIds = c.Ints(nf);
        two.mHasIds = true;
    }
    if (flags[4]) {
        for (int k = 0; k < 3; ++k)
            two.mHierarchy.push_back(c.Ints(np2));
        two.mHasHierarchy = true;
    }

    AnimFamily three, one, sph;
    bool has3 = false, has1 = false, has_sph = false;
    if (flags[2]) {
        const std::vector<std::int64_t> h = c.Ints(4);  // elements, parts, scalars, tensors
        three.mKey = "3d";
        three.mNodesPer = 8;
        three.mConn = c.Ints(8 * h[0]);
        anim_section(c, flags, three, h[1], h[2], 6, h[3], false);
        has3 = true;
    }
    if (flags[3]) {
        const std::vector<std::int64_t> h = c.Ints(5);  // elements, parts, scalars, torsors, skews
        one.mKey = "1d";
        one.mNodesPer = 2;
        one.mConn = c.Ints(2 * h[0]);
        anim_section(c, flags, one, h[1], h[2], 9, h[3], h[4] != 0);
        has1 = true;
    }
    if (flags[4]) {
        const std::int64_t nsub = c.Int();
        for (std::int64_t k = 0; k < nsub; ++k) {
            c.Text(50);
            c.Int();          // parent
            c.Ints(c.Int());  // child subsets
            for (int d = 0; d < 3; ++d)
                c.Ints(c.Int());  // 2-D, 3-D and 1-D parts
        }
        const std::int64_t nmat = c.Int(), nprop = c.Int();
        c.Texts(nmat, 50);
        c.Ints(nmat);
        c.Texts(nprop, 50);
        c.Ints(nprop);
    }
    if (flags[5]) {
        const std::vector<std::int64_t> th = c.Ints(4);  // nodes, 2-D, 3-D, 1-D elements
        for (std::int64_t count : th) {
            c.Ints(count);
            c.Texts(count, 50);
        }
    }
    if (flags[7]) {
        const std::vector<std::int64_t> h = c.Ints(4);  // elements, parts, scalars, tensors
        sph.mKey = "sph";
        sph.mNodesPer = 1;
        if (h[0]) {
            sph.mConn = c.Ints(h[0]);
            sph.mAlive = c.Flags(h[0]);
        }
        if (h[1]) {
            sph.mPartEnds = c.Ints(h[1]);
            sph.mPartTexts = c.Texts(h[1], 50);
        }
        if (h[2]) {
            const std::vector<std::string> names = c.Texts(h[2], 81);
            const std::vector<double> values = c.Floats(h[2] * h[0]);
            for (std::int64_t k = 0; k < h[2]; ++k)
                sph.mScalars.emplace_back(names[static_cast<std::size_t>(k)],
                                          std::vector<double>(values.begin() + k * h[0],
                                                              values.begin() + (k + 1) * h[0]));
        }
        if (h[3]) {
            const std::vector<std::string> names = c.Texts(h[3], 81);
            const std::vector<double> values = c.Floats(h[0] * h[3] * 6);
            for (std::int64_t k = 0; k < h[3]; ++k)
                sph.mTensors.emplace_back(names[static_cast<std::size_t>(k)], 6,
                                          std::vector<double>(values.begin() + k * h[0] * 6,
                                                              values.begin() + (k + 1) * h[0] * 6));
        }
        if (flags[0] == 1) {
            sph.mMass = c.Floats(h[0]);
            sph.mHasMass = true;
        }
        if (flags[1] == 1) {
            sph.mIds = c.Ints(h[0]);
            sph.mHasIds = true;
        }
        if (flags[4]) {
            for (int k = 0; k < 3; ++k)
                sph.mHierarchy.push_back(c.Ints(h[1]));
            sph.mHasHierarchy = true;
        }
        has_sph = true;
    }
    // The families in file order: 1-D, 2-D, 3-D, SPH.
    if (has1 && one.Size())
        families.push_back(std::move(one));
    if (two.Size())
        families.push_back(std::move(two));
    if (has3 && three.Size())
        families.push_back(std::move(three));
    if (has_sph && sph.Size())
        families.push_back(std::move(sph));

    // Cells: each element's type and nodes, then one block per type.
    struct AnimCell {
        std::size_t mFamily, mElement;
        std::vector<std::int64_t> mNodes;
    };
    std::vector<std::string> block_types;
    std::map<std::string, std::vector<AnimCell>> members;
    std::size_t skipped = 0;
    for (std::size_t fi = 0; fi < families.size(); ++fi) {
        const AnimFamily& fam = families[fi];
        for (std::size_t e = 0; e < fam.Size(); ++e) {
            std::vector<std::int64_t> row(
                fam.mConn.begin() + static_cast<std::ptrdiff_t>(e * fam.mNodesPer),
                fam.mConn.begin() + static_cast<std::ptrdiff_t>((e + 1) * fam.mNodesPer));
            for (std::int64_t v : row)
                if (v < 0 || v >= nn)
                    throw ReadError("Radioss animation: '" + rPath + "' names a node out of range");
            std::string type;
            if (fam.mKey == "1d") {
                type = "line";
            } else if (fam.mKey == "sph") {
                type = "vertex";
            } else if (fam.mKey == "2d") {
                std::vector<std::int64_t> distinct;
                for (std::int64_t v : row)
                    if (std::find(distinct.begin(), distinct.end(), v) == distinct.end())
                        distinct.push_back(v);
                row = distinct;
                type = row.size() == 4 ? "quad" : (row.size() == 3 ? "triangle" : "");
            } else {
                std::array<std::int64_t, 8> a{};
                std::copy(row.begin(), row.end(), a.begin());
                detail::CollapsedBrick b = detail::collapse_brick(a);
                type = b.mType;
                row = std::move(b.mNodes);
            }
            if (type.empty()) {
                ++skipped;
                continue;
            }
            auto [it, fresh] = members.emplace(type, std::vector<AnimCell>{});
            if (fresh)
                block_types.push_back(type);
            it->second.push_back({fi, e, std::move(row)});
        }
    }
    if (skipped)
        log::warn("Radioss animation: {} degenerate facet(s) skipped", skipped);

    Mesh mesh;
    NDArray points(DType::Float64, {static_cast<std::size_t>(nn), 3});
    std::copy(coords.begin(), coords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    std::map<std::pair<std::size_t, std::size_t>, std::int64_t> cell_of;
    std::vector<int> cell_dim;
    for (const std::string& type : block_types) {
        const std::vector<AnimCell>& list = members[type];
        const std::size_t k = list.front().mNodes.size();
        NDArray conn(DType::Int64, {list.size(), k});
        for (std::size_t r = 0; r < list.size(); ++r) {
            std::copy(list[r].mNodes.begin(), list[r].mNodes.end(),
                      conn.As<std::int64_t>() + r * k);
            cell_of[{list[r].mFamily, list[r].mElement}] =
                static_cast<std::int64_t>(cell_dim.size());
            cell_dim.push_back(type == "vertex"                         ? 0
                               : type == "line"                         ? 1
                               : (type == "triangle" || type == "quad") ? 2
                                                                        : 3);
        }
        mesh.AddCellBlock(type, std::move(conn));
    }
    mesh.AddFieldData("meshio:time", [&] {
        NDArray a(DType::Float64, {1});
        a.As<double>()[0] = time;
        return a;
    }());

    auto vec = [](const std::vector<double>& rV, std::vector<std::size_t> Shape) {
        NDArray a(DType::Float64, std::move(Shape));
        std::copy(rV.begin(), rV.end(), a.As<double>());
        return a;
    };
    const std::size_t n_points = static_cast<std::size_t>(nn);
    for (const auto& [name, values] : point_scalars)
        mesh.AddPointData(name, vec(values, {n_points}));
    for (std::size_t k = 0; k < vnames.size(); ++k)
        mesh.AddPointData(
            vnames[k],
            vec(std::vector<double>(
                    vvalues.begin() + static_cast<std::ptrdiff_t>(k * 3 * n_points),
                    vvalues.begin() + static_cast<std::ptrdiff_t>((k + 1) * 3 * n_points)),
                {n_points, 3}));
    if (flags[1]) {
        NDArray a(DType::Int64, {n_points});
        std::copy(node_ids.begin(), node_ids.end(), a.As<std::int64_t>());
        mesh.AddPointData("radioss:node_id", std::move(a));
    }
    if (flags[0] == 1)
        mesh.AddPointData("radioss:mass", vec(point_mass, {n_points}));
    if (block_types.empty())
        return mesh;

    // One array per block from `rGet(family, element, out)`.
    auto per_block = [&](std::size_t Width, double Fill, DType Type, auto&& rGet) {
        std::vector<NDArray> out;
        for (const std::string& type : block_types) {
            const std::vector<AnimCell>& list = members[type];
            const std::size_t w = Width ? Width : 1;
            NDArray a(Type, Width ? std::vector<std::size_t>{list.size(), Width}
                                  : std::vector<std::size_t>{list.size()});
            std::vector<double> row(w);
            for (std::size_t r = 0; r < list.size(); ++r) {
                std::fill(row.begin(), row.end(), Fill);
                rGet(list[r].mFamily, list[r].mElement, row.data());
                for (std::size_t j = 0; j < w; ++j) {
                    if (Type == DType::Float64)
                        a.As<double>()[r * w + j] = row[j];
                    else if (Type == DType::Int8)
                        a.As<std::int8_t>()[r * w + j] = static_cast<std::int8_t>(row[j]);
                    else
                        a.As<std::int64_t>()[r * w + j] = static_cast<std::int64_t>(row[j]);
                }
            }
            out.push_back(std::move(a));
        }
        return out;
    };

    // Per family and element: the part id, and the part's hierarchy values.
    std::vector<std::vector<std::int64_t>> part_id(families.size());
    std::vector<std::vector<std::array<std::int64_t, 2>>> part_matprop(families.size());
    struct AnimPart {
        std::string mName;
        std::int64_t mId;
        std::size_t mFamily;
        std::int64_t mStart, mEnd;
    };
    std::vector<AnimPart> parts;
    for (std::size_t fi = 0; fi < families.size(); ++fi) {
        const AnimFamily& fam = families[fi];
        part_id[fi].assign(fam.Size(), 0);
        part_matprop[fi].assign(fam.Size(), {0, 0});
        anim_for_parts(fam, [&](std::size_t k, std::int64_t b, std::int64_t e) {
            const std::string& text = fam.mPartTexts[k];
            const std::size_t colon = text.find(':');
            std::int64_t pid = 0;
            try {
                pid = std::stoll(text.substr(0, colon));
            } catch (...) {
                pid = 0;
            }
            std::string name = colon == std::string::npos ? text : text.substr(colon + 1);
            const std::size_t nb = name.find_first_not_of(" \t");
            name = nb == std::string::npos ? std::string()
                                           : name.substr(nb, name.find_last_not_of(" \t") - nb + 1);
            for (std::int64_t i = b; i < e; ++i) {
                part_id[fi][static_cast<std::size_t>(i)] = pid;
                if (fam.mHasHierarchy)
                    part_matprop[fi][static_cast<std::size_t>(i)] = {fam.mHierarchy[1][k],
                                                                     fam.mHierarchy[2][k]};
            }
            parts.push_back({name, pid, fi, b, e});
        });
    }
    mesh.AddCellData("radioss:part",
                     per_block(0, 0.0, DType::Int64, [&](std::size_t f, std::size_t e, double* p) {
                         *p = static_cast<double>(part_id[f][e]);
                     }));
    // A deleted (eroded) element's flag is 0; any other byte is alive (the
    // facets write 0xFF, bricks and beams 1).
    mesh.AddCellData("radioss:alive",
                     per_block(0, 0.0, DType::Int8, [&](std::size_t f, std::size_t e, double* p) {
                         *p = families[f].mAlive[e] ? 1.0 : 0.0;
                     }));
    auto all = [&](bool AnimFamily::* pMember) {
        for (const AnimFamily& fam : families)
            if (!(fam.*pMember))
                return false;
        return true;
    };
    if (all(&AnimFamily::mHasIds))
        mesh.AddCellData(
            "radioss:element_id",
            per_block(0, 0.0, DType::Int64, [&](std::size_t f, std::size_t e, double* p) {
                *p = static_cast<double>(families[f].mIds[e]);
            }));
    if (all(&AnimFamily::mHasMass))
        mesh.AddCellData("radioss:mass", per_block(0, 0.0, DType::Float64,
                                                   [&](std::size_t f, std::size_t e, double* p) {
                                                       *p = families[f].mMass[e];
                                                   }));
    if (all(&AnimFamily::mHasHierarchy)) {
        mesh.AddCellData(
            "radioss:material",
            per_block(0, 0.0, DType::Int64, [&](std::size_t f, std::size_t e, double* p) {
                *p = static_cast<double>(part_matprop[f][e][0]);
            }));
        mesh.AddCellData(
            "radioss:property",
            per_block(0, 0.0, DType::Int64, [&](std::size_t f, std::size_t e, double* p) {
                *p = static_cast<double>(part_matprop[f][e][1]);
            }));
    }

    // Element results: the union of the families' names, NaN where a family has
    // none; tensors as xx yy zz xy yz zx (a 2-D one's out-of-plane parts NaN),
    // 1-D force/moment sets with their nine components.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, std::map<std::size_t, const std::vector<double>*>> scalars;
    std::map<std::string, std::map<std::size_t, std::pair<std::size_t, const std::vector<double>*>>>
        tensors;
    for (std::size_t fi = 0; fi < families.size(); ++fi) {
        for (const auto& [name, values] : families[fi].mScalars)
            scalars[name][fi] = &values;
        for (const auto& [name, width, values] : families[fi].mTensors)
            tensors[name][fi] = {width, &values};
    }
    for (const auto& [name, by_family] : scalars)
        mesh.AddCellData(
            name, per_block(0, nan, DType::Float64, [&](std::size_t f, std::size_t e, double* p) {
                const auto it = by_family.find(f);
                if (it != by_family.end())
                    *p = (*it->second)[e];
            }));
    for (const auto& [name, by_family] : tensors) {
        bool nine = false, other = false;
        for (const auto& [f, wv] : by_family)
            (wv.first == 9 ? nine : other) = true;
        if (nine && other)
            throw ReadError("Radioss animation: '" + rPath +
                            "' names both a 1-D force set and a tensor '" + name + "'");
        const std::size_t width = nine ? 9 : 6;
        mesh.AddCellData(name, per_block(width, nan, DType::Float64,
                                         [&](std::size_t f, std::size_t e, double* p) {
                                             const auto it = by_family.find(f);
                                             if (it == by_family.end())
                                                 return;
                                             const std::size_t w = it->second.first;
                                             const double* v = it->second.second->data() + e * w;
                                             if (w == 3) {  // xx yy xy
                                                 p[0] = v[0];
                                                 p[1] = v[1];
                                                 p[3] = v[2];
                                             } else {
                                                 std::copy(v, v + w, p);
                                             }
                                         }));
    }

    // Parts -> cell regions (a part id can name a 1-D and a 2-D part).
    std::map<std::pair<std::string, std::int64_t>, std::vector<std::int64_t>> regions;
    std::vector<std::pair<std::string, std::int64_t>> order;
    for (const AnimPart& p : parts) {
        const std::pair<std::string, std::int64_t> key{
            p.mName.empty() ? "Part " + std::to_string(p.mId) : p.mName, p.mId};
        auto [it, fresh] = regions.emplace(key, std::vector<std::int64_t>{});
        if (fresh)
            order.push_back(key);
        for (std::int64_t i = p.mStart; i < p.mEnd; ++i) {
            const auto c = cell_of.find({p.mFamily, static_cast<std::size_t>(i)});
            if (c != cell_of.end())
                it->second.push_back(c->second);
        }
    }
    for (const auto& key : order) {
        std::vector<std::int64_t>& ids = regions[key];
        std::sort(ids.begin(), ids.end());
        int dim = -1;
        for (std::int64_t c : ids)
            dim = std::max(dim, cell_dim[static_cast<std::size_t>(c)]);
        NDArray a(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), a.As<std::int64_t>());
        mesh.AddRegion(Region(key.first, RegionKind::Cell, dim, key.second, std::move(a)));
    }
    return mesh;
}

MeshMetadata read_radioss_anim_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    const Mesh mesh = read_radioss_anim(rPath);
    MeshMetadata meta = metadata_from_mesh(mesh);
    meta.mFellBackToFullRead = true;
    meta.mFormat = "radioss_anim";
    meta.mTimeValues = {mesh.FieldData("meshio:time").As<double>()[0]};
    return meta;
}

}  // namespace meshioplusplus
