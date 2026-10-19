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
#include <memory>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/radioss_anim.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/degenerate_solid.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

constexpr std::int32_t kAnimMagic = 0x542C;
/// Region names of the time-history items: the file's label behind this prefix.
constexpr const char* kAnimThPrefix = "radioss_th:";

// Big-endian reads with a bounds check.
class AnimCursor {
public:
    AnimCursor(std::string_view rData, const std::string& rPath) : mData(rData), mPath(rPath) {}

    std::size_t Remaining() const { return mData.size() - mPos; }

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
        if (Count > 0 && N > 0)  // bounded by the bytes left, as Take() checks
            out.reserve(std::min<std::size_t>(static_cast<std::size_t>(Count), Remaining() / N));
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
    // Time-history items: per category (nodes, 2-D, 3-D, 1-D elements) the 0-based index of a
    // node or of an element of that family, each with a label.
    std::array<std::vector<std::pair<std::int64_t, std::string>>, 4> th_items;
    if (flags[5]) {
        const std::vector<std::int64_t> th = c.Ints(4);
        for (std::size_t slot = 0; slot < 4; ++slot) {
            const std::vector<std::int64_t> ids = c.Ints(th[slot]);
            const std::vector<std::string> texts = c.Texts(th[slot], 50);
            for (std::size_t k = 0; k < ids.size(); ++k)
                th_items[slot].emplace_back(ids[k], texts[k]);
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
    // The time-history items become regions named after their label: a point region for the
    // nodes, cell regions for the elements (a label can name elements of several families).
    const auto add_th_regions = [&]() {
        const auto label = [](const std::string& rText) {
            const std::size_t b = rText.find_first_not_of(" \t");
            return std::string(kAnimThPrefix) +
                   (b == std::string::npos
                        ? std::string()
                        : rText.substr(b, rText.find_last_not_of(" \t") - b + 1));
        };
        std::map<std::string, std::set<std::int64_t>> nodes, cells;
        std::vector<std::string> node_order, cell_order;
        for (const auto& [index, text] : th_items[0]) {
            if (index < 0 || index >= nn)
                continue;
            const std::string name = label(text);
            if (nodes.try_emplace(name).second)
                node_order.push_back(name);
            nodes[name].insert(index);
        }
        static constexpr const char* kKeys[3] = {"2d", "3d", "1d"};
        for (std::size_t slot = 1; slot < 4; ++slot) {
            for (std::size_t fi = 0; fi < families.size(); ++fi) {
                if (families[fi].mKey != kKeys[slot - 1])
                    continue;
                for (const auto& [index, text] : th_items[slot]) {
                    const auto c = cell_of.find({fi, static_cast<std::size_t>(index)});
                    if (index < 0 || c == cell_of.end())
                        continue;
                    const std::string name = label(text);
                    if (cells.try_emplace(name).second)
                        cell_order.push_back(name);
                    cells[name].insert(c->second);
                }
            }
        }
        const auto region = [&](const std::string& rName, RegionKind Kind,
                                const std::set<std::int64_t>& rIds, int Dim) {
            NDArray a(DType::Int64, {rIds.size()});
            std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
            mesh.AddRegion(Region(rName, Kind, Dim, -1, std::move(a)));
        };
        for (const std::string& name : node_order)
            region(name, RegionKind::Point, nodes[name], -1);
        for (const std::string& name : cell_order) {
            int dim = -1;
            for (const std::int64_t c : cells[name])
                dim = std::max(dim, cell_dim[static_cast<std::size_t>(c)]);
            region(name, RegionKind::Cell, cells[name], dim);
        }
    };
    if (block_types.empty()) {
        add_th_regions();
        return mesh;
    }

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
    add_th_regions();
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


// ---------------------------------------------------------------------------
// The writer: the reader's inverse. One mesh is one state; a series is one file
// per step. The Python twin (`radioss_anim/_write.py`) writes the same bytes.
// ---------------------------------------------------------------------------
namespace {

/// Big-endian output buffer.
class AwOut {
public:
    void Ints(const std::vector<std::int64_t>& rValues) {
        for (const std::int64_t v : rValues) {
            const std::uint32_t u = static_cast<std::uint32_t>(static_cast<std::int32_t>(v));
            Be32(u);
        }
    }
    void Floats(const std::vector<double>& rValues) {
        for (const double v : rValues)
            Float(v);
    }
    void Float(double Value) {
        const float f = static_cast<float>(Value);
        std::uint32_t u;
        std::memcpy(&u, &f, 4);
        Be32(u);
    }
    void Be32(std::uint32_t U) {
        mBytes.push_back(static_cast<char>(U >> 24));
        mBytes.push_back(static_cast<char>(U >> 16));
        mBytes.push_back(static_cast<char>(U >> 8));
        mBytes.push_back(static_cast<char>(U));
    }
    void Text(const std::string& rText, std::size_t Width) {
        const std::string cut = rText.substr(0, Width - 1);
        mBytes += cut;
        mBytes.append(Width - cut.size(), '\0');
    }
    void Bytes(const std::vector<std::uint8_t>& rValues) {
        mBytes.append(reinterpret_cast<const char*>(rValues.data()), rValues.size());
    }
    void Zeros(std::size_t N) { mBytes.append(N, '\0'); }
    void Append(const AwOut& rOther) { mBytes += rOther.mBytes; }
    const std::string& Data() const { return mBytes; }

private:
    std::string mBytes;
};

struct AwFamily {
    std::vector<std::size_t> mIndex;                // global cell numbers, ordered by part
    std::vector<std::vector<std::int64_t>> mRows;  // file node rows
    std::vector<std::int64_t> mParts;               // part id per element
    std::size_t Size() const { return mIndex.size(); }
};

/// A cell array over all blocks (block-major) as doubles, with its width.
struct AwColumn {
    std::string mName;
    std::vector<double> mValues;
    std::size_t mWidth = 1;
};

enum AwKey : std::size_t { kAw2d = 0, kAw3d = 1, kAw1d = 2, kAwSph = 3 };

std::string aw_pad_id(std::int64_t Id) {
    std::string s = std::to_string(Id);
    if (s.size() < 9)
        s.insert(0, 9 - s.size(), ' ');
    return s;
}

std::vector<double> aw_values(const NDArray& rArray) {
    const detail::DoubleView view(rArray);
    return std::vector<double>(view.Data(), view.Data() + rArray.Size());
}

std::vector<std::int64_t> aw_int_values(const NDArray& rArray) {
    const detail::Int64View view(rArray);
    return std::vector<std::int64_t>(view.Data(), view.Data() + rArray.Size());
}

/// The block-major cell array `rName` flattened, or false.
bool aw_concat(const Mesh& rMesh, const std::string& rName, std::vector<double>& rOut) {
    if (!rMesh.HasCellData(rName))
        return false;
    rOut.clear();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::vector<double> v = aw_values(rMesh.CellData(rName, b));
        rOut.insert(rOut.end(), v.begin(), v.end());
    }
    return true;
}

bool aw_concat_int(const Mesh& rMesh, const std::string& rName, std::vector<std::int64_t>& rOut) {
    if (!rMesh.HasCellData(rName))
        return false;
    rOut.clear();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::vector<std::int64_t> v = aw_int_values(rMesh.CellData(rName, b));
        rOut.insert(rOut.end(), v.begin(), v.end());
    }
    return true;
}

int aw_family_of(const std::string& rType) {
    if (rType == "vertex")
        return kAwSph;
    if (rType == "line")
        return kAw1d;
    if (rType == "triangle" || rType == "quad")
        return kAw2d;
    if (rType == "tetra" || rType == "pyramid" || rType == "wedge" || rType == "hexahedron")
        return kAw3d;
    return -1;
}

bool aw_special(const std::string& rName) {
    return rName == "radioss:part" || rName == "radioss:alive" || rName == "radioss:element_id" ||
           rName == "radioss:mass" || rName == "radioss:material" || rName == "radioss:property";
}

double aw_time_of(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kSequenceTimeKey))
        return 0.0;
    const NDArray& a = rMesh.FieldData(kSequenceTimeKey);
    return a.Size() == 0 ? 0.0 : detail::DoubleView(a).Data()[0];
}

std::string aw_state_bytes(const Mesh& rMesh, double Time) {
    const std::size_t nn = rMesh.NumPoints();
    std::vector<double> points(nn * 3, 0.0);
    {
        const NDArray& pts = rMesh.Points();
        const detail::DoubleView view(pts);
        const std::size_t width = nn == 0 ? 3 : pts.Size() / nn;
        for (std::size_t i = 0; i < nn; ++i)
            for (std::size_t k = 0; k < 3 && k < width; ++k)
                points[i * 3 + k] = view.Data()[i * width + k];
    }
    // The element families.
    std::array<AwFamily, 4> fam;
    std::map<std::string, std::size_t> skipped;
    std::size_t ncells = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t count = cb.NumCells();
        const int key = aw_family_of(cb.Type());
        if (key < 0) {
            if (count > 0)
                skipped[cb.Type()] += count;
        } else {
            const detail::Int64View conn(cb.Conn());
            const std::size_t w = cb.NodesPerCell();
            AwFamily& f = fam[static_cast<std::size_t>(key)];
            for (std::size_t r = 0; r < count; ++r) {
                const std::int64_t* row = conn.Data() + r * w;
                for (std::size_t k = 0; k < w; ++k)
                    if (row[k] < 0 || static_cast<std::size_t>(row[k]) >= nn)
                        throw WriteError("Radioss animation: a '" + cb.Type() +
                                         "' cell refers to a point the mesh does not have");
                f.mIndex.push_back(ncells + r);
                if (key == kAwSph) {
                    f.mRows.push_back({row[0]});
                } else if (key == kAw2d) {
                    std::vector<std::int64_t> v(row, row + w);
                    if (w == 3)
                        v.push_back(row[2]);
                    f.mRows.push_back(std::move(v));
                } else if (key == kAw3d) {
                    const auto brick = detail::expand_brick(cb.Type(), row);
                    f.mRows.emplace_back(brick.begin(), brick.end());
                } else {
                    f.mRows.push_back({row[0], row[1]});
                }
            }
        }
        ncells += count;
    }
    for (const auto& [type, count] : skipped)
        log::warn("Radioss animation: {} '{}' cell(s) are not written", count, type);

    // Per-cell attributes.
    std::vector<std::int64_t> part_arr, eids, emat, eprop;
    std::vector<double> alive_d, emass;
    const bool have_part = aw_concat_int(rMesh, "radioss:part", part_arr);
    const bool have_alive = aw_concat(rMesh, "radioss:alive", alive_d);
    const bool have_eids = aw_concat_int(rMesh, "radioss:element_id", eids);
    const bool have_emass = aw_concat(rMesh, "radioss:mass", emass);
    const bool have_emat = aw_concat_int(rMesh, "radioss:material", emat);
    const bool have_eprop = aw_concat_int(rMesh, "radioss:property", eprop);
    std::vector<std::int64_t> node_ids;
    std::vector<double> nmass;
    const bool have_nids = rMesh.HasPointData("radioss:node_id");
    const bool have_nmass = rMesh.HasPointData("radioss:mass");
    if (have_nids)
        node_ids = aw_int_values(rMesh.PointData("radioss:node_id"));
    if (have_nmass)
        nmass = aw_values(rMesh.PointData("radioss:mass"));

    // Parts: radioss:part, else the cell regions, else part 1.
    std::vector<std::int64_t> region_pid(ncells, -1);
    std::map<std::int64_t, std::string> names_by_pid;
    bool have_regions = false;
    {
        std::size_t k = 0;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const auto& reg = rMesh.Region(i);
            if (reg.mKind != RegionKind::Cell ||
                reg.mName.rfind(kAnimThPrefix, 0) == 0)  // time-history items are not parts
                continue;
            have_regions = true;
            const std::int64_t pid = reg.mTag >= 0 ? reg.mTag : static_cast<std::int64_t>(k + 1);
            ++k;
            names_by_pid.emplace(pid, reg.mName);
            if (!have_part) {
                const detail::Int64View entries(reg.mEntries);
                for (std::size_t e = 0; e < reg.mEntries.Size(); ++e) {
                    const std::int64_t c = entries.Data()[e];
                    if (c >= 0 && static_cast<std::size_t>(c) < ncells &&
                        region_pid[static_cast<std::size_t>(c)] < 0)
                        region_pid[static_cast<std::size_t>(c)] = pid;
                }
            }
        }
    }
    auto pid_of = [&](std::size_t c) -> std::int64_t {
        if (have_part)
            return part_arr[c];
        if (have_regions && region_pid[c] >= 0)
            return region_pid[c];
        return 1;
    };
    for (AwFamily& f : fam) {
        std::vector<std::size_t> order(f.Size());
        for (std::size_t k = 0; k < order.size(); ++k)
            order[k] = k;
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return pid_of(f.mIndex[a]) < pid_of(f.mIndex[b]);
        });
        AwFamily sorted;
        for (const std::size_t k : order) {
            sorted.mIndex.push_back(f.mIndex[k]);
            sorted.mRows.push_back(std::move(f.mRows[k]));
            sorted.mParts.push_back(pid_of(f.mIndex[k]));
        }
        f = std::move(sorted);
    }

    // Time-history items: each entry of a `radioss_th:<label>` region, as the index of a node or
    // of an element in its family (the order the sections above were just sorted into).
    std::array<std::vector<std::pair<std::int64_t, std::string>>, 4> th_items;
    {
        std::map<std::size_t, std::pair<std::size_t, std::size_t>> position;  // cell -> family, row
        for (std::size_t key = 0; key < fam.size(); ++key)
            for (std::size_t k = 0; k < fam[key].Size(); ++k)
                position[fam[key].mIndex[k]] = {key, k};
        const std::string prefix(kAnimThPrefix);
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const auto& reg = rMesh.Region(i);
            if (reg.mName.rfind(prefix, 0) != 0)
                continue;
            const std::string label = reg.mName.substr(prefix.size());
            const detail::Int64View entries(reg.mEntries);
            for (std::size_t e = 0; e < reg.mEntries.Size(); ++e) {
                const std::int64_t v = entries.Data()[e];
                if (reg.mKind == RegionKind::Point) {
                    if (v >= 0 && static_cast<std::size_t>(v) < nn)
                        th_items[0].emplace_back(v, label);
                } else if (reg.mKind == RegionKind::Cell) {
                    const auto at = position.find(static_cast<std::size_t>(v));
                    if (v < 0 || at == position.end() || at->second.first == kAwSph)
                        continue;
                    const std::size_t slot = at->second.first == kAw2d   ? 1
                                             : at->second.first == kAw3d ? 2
                                                                         : 3;
                    th_items[slot].emplace_back(static_cast<std::int64_t>(at->second.second),
                                                label);
                }
            }
        }
    }
    bool have_th = false;
    for (const auto& items : th_items)
        have_th = have_th || !items.empty();

    // Nodal and element results.
    std::vector<AwColumn> scalars, vectors, escalars, tensors, forces;
    for (const std::string& name : rMesh.PointDataNames()) {
        if (name == "radioss:node_id" || name == "radioss:mass")
            continue;
        AwColumn col;
        col.mName = name;
        col.mValues = aw_values(rMesh.PointData(name));
        if (nn > 0 && col.mValues.size() == nn)
            scalars.push_back(std::move(col));
        else if (nn > 0 && col.mValues.size() == 3 * nn)
            vectors.push_back(std::move(col));
        else
            log::warn("Radioss animation: point data '{}' is not written", name);
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        if (aw_special(name))
            continue;
        AwColumn col;
        col.mName = name;
        aw_concat(rMesh, name, col.mValues);
        if (ncells == 0 || col.mValues.size() % ncells != 0) {
            log::warn("Radioss animation: cell data '{}' is not written", name);
            continue;
        }
        col.mWidth = col.mValues.size() / ncells;
        if (col.mWidth == 1)
            escalars.push_back(std::move(col));
        else if (col.mWidth == 6)
            tensors.push_back(std::move(col));
        else if (col.mWidth == 9)
            forces.push_back(std::move(col));
        else
            log::warn("Radioss animation: cell data '{}' is not written", name);
    }

    const bool have_mass = have_emass || have_nmass;
    const bool have_ids = have_eids || have_nids;
    const bool have_hier = have_emat || have_eprop;
    const std::vector<std::int64_t> flags = {have_mass ? 1 : 0,
                                             have_ids ? 1 : 0,
                                             fam[kAw3d].Size() ? 1 : 0,
                                             fam[kAw1d].Size() ? 1 : 0,
                                             have_hier ? 1 : 0,
                                             have_th ? 1 : 0,
                                             0,
                                             fam[kAwSph].Size() ? 1 : 0,
                                             0,
                                             0};

    struct PartTable {
        std::vector<std::int64_t> mEnds, mMats, mProps;
        std::vector<std::string> mTexts;
    };
    auto part_table = [&](const AwFamily& f) {
        PartTable t;
        for (std::size_t k = 0; k < f.mParts.size(); ++k) {
            if (k + 1 != f.mParts.size() && f.mParts[k + 1] == f.mParts[k])
                continue;
            const std::int64_t pid = f.mParts[k];
            t.mEnds.push_back(static_cast<std::int64_t>(k + 1));
            const auto name = names_by_pid.find(pid);
            t.mTexts.push_back(aw_pad_id(pid) + ":" +
                               (name != names_by_pid.end() ? name->second
                                                           : "Part " + std::to_string(pid)));
            const std::size_t c = f.mIndex[k];
            t.mMats.push_back(have_emat ? emat[c] : 0);
            t.mProps.push_back(have_eprop ? eprop[c] : 0);
        }
        return t;
    };
    auto alive_bytes = [&](const AwFamily& f, bool Facets) {
        std::vector<std::uint8_t> a;
        for (const std::size_t c : f.mIndex) {
            const bool on = have_alive ? alive_d[c] != 0.0 : true;
            a.push_back(on ? (Facets ? 0xFF : 1) : 0);
        }
        return a;
    };
    auto masses_of = [&](const AwFamily& f) {
        std::vector<double> m;
        for (const std::size_t c : f.mIndex)
            m.push_back(have_emass ? emass[c] : 0.0);
        return m;
    };
    auto ids_of = [&](const AwFamily& f) {
        std::vector<std::int64_t> v;
        for (const std::size_t c : f.mIndex)
            v.push_back(have_eids ? eids[c] : static_cast<std::int64_t>(c + 1));
        return v;
    };
    auto rows_flat = [](const AwFamily& f) {
        std::vector<std::int64_t> v;
        for (const auto& r : f.mRows)
            v.insert(v.end(), r.begin(), r.end());
        return v;
    };
    // `width` columns of `rCol` for each element of a family.
    auto gather = [&](const AwColumn& rCol, const AwFamily& f, const std::vector<int>& rCols) {
        std::vector<double> v;
        for (const std::size_t c : f.mIndex)
            for (const int k : rCols)
                v.push_back(rCol.mValues[c * rCol.mWidth + static_cast<std::size_t>(k)]);
        return v;
    };

    AwOut out;
    out.Ints({kAnimMagic});
    out.Float(Time);
    for (const char* t : {"Time=", "ModAnim", "Radioss Run="})
        out.Text(t, 81);
    out.Ints(flags);

    // 2-D: the nodes and the facets.
    const AwFamily& f2 = fam[kAw2d];
    const PartTable pt2 = part_table(f2);
    out.Ints({static_cast<std::int64_t>(nn), static_cast<std::int64_t>(f2.Size()),
              static_cast<std::int64_t>(pt2.mEnds.size()),
              static_cast<std::int64_t>(scalars.size()),
              static_cast<std::int64_t>(escalars.size()),
              static_cast<std::int64_t>(vectors.size()),
              static_cast<std::int64_t>(tensors.size()), 0});
    out.Floats(points);
    out.Ints(rows_flat(f2));
    out.Bytes(alive_bytes(f2, true));
    out.Ints(pt2.mEnds);
    for (const std::string& t : pt2.mTexts)
        out.Text(t, 50);
    out.Zeros(2 * 3 * nn);  // nodal normals
    if (!scalars.empty() || !escalars.empty()) {
        for (const AwColumn& c : scalars)
            out.Text(c.mName, 81);
        for (const AwColumn& c : escalars)
            out.Text(c.mName, 81);
        for (const AwColumn& c : scalars)
            out.Floats(c.mValues);
        for (const AwColumn& c : escalars)
            out.Floats(gather(c, f2, {0}));
    }
    for (const AwColumn& c : vectors)
        out.Text(c.mName, 81);
    for (const AwColumn& c : vectors)
        out.Floats(c.mValues);
    if (!tensors.empty()) {
        for (const AwColumn& c : tensors)
            out.Text(c.mName, 81);
        for (const AwColumn& c : tensors)
            out.Floats(gather(c, f2, {0, 1, 3}));
    }
    if (flags[0]) {
        out.Floats(masses_of(f2));
        std::vector<double> pm(nn, 0.0);
        if (have_nmass)
            for (std::size_t i = 0; i < nn && i < nmass.size(); ++i)
                pm[i] = nmass[i];
        out.Floats(pm);
    }
    if (flags[1]) {
        std::vector<std::int64_t> nid(nn);
        for (std::size_t i = 0; i < nn; ++i)
            nid[i] = have_nids ? node_ids[i] : static_cast<std::int64_t>(i + 1);
        out.Ints(nid);
        out.Ints(ids_of(f2));
    }
    if (flags[4]) {
        out.Ints(std::vector<std::int64_t>(pt2.mEnds.size(), 1));
        out.Ints(pt2.mMats);
        out.Ints(pt2.mProps);
    }

    std::array<std::size_t, 4> part_counts = {pt2.mEnds.size(), 0, 0, 0};
    // The 3-D, 1-D and SPH sections share a tail.
    auto section = [&](const AwFamily& f, const std::vector<std::int64_t>& rCounts,
                       const std::vector<AwColumn>& rScalars,
                       const std::vector<AwColumn>& rTensors, bool Facets, std::size_t Slot) {
        const PartTable pt = part_table(f);
        part_counts[Slot] = pt.mEnds.size();
        std::vector<std::int64_t> counts = rCounts;
        counts[1] = static_cast<std::int64_t>(pt.mEnds.size());
        out.Ints(counts);
        out.Ints(rows_flat(f));
        out.Bytes(alive_bytes(f, Facets));
        out.Ints(pt.mEnds);
        for (const std::string& t : pt.mTexts)
            out.Text(t, 50);
        if (!rScalars.empty()) {
            for (const AwColumn& c : rScalars)
                out.Text(c.mName, 81);
            for (const AwColumn& c : rScalars)
                out.Floats(gather(c, f, {0}));
        }
        if (!rTensors.empty()) {
            for (const AwColumn& c : rTensors)
                out.Text(c.mName, 81);
            std::vector<int> all(rTensors.front().mWidth);
            for (std::size_t k = 0; k < all.size(); ++k)
                all[k] = static_cast<int>(k);
            for (const AwColumn& c : rTensors)
                out.Floats(gather(c, f, all));
        }
        if (flags[0])
            out.Floats(masses_of(f));
        if (flags[1])
            out.Ints(ids_of(f));
        if (flags[4]) {
            out.Ints(std::vector<std::int64_t>(pt.mEnds.size(), 1));
            out.Ints(pt.mMats);
            out.Ints(pt.mProps);
        }
    };
    const AwFamily& f3 = fam[kAw3d];
    if (flags[2])
        section(f3,
                {static_cast<std::int64_t>(f3.Size()), 0,
                 static_cast<std::int64_t>(escalars.size()),
                 static_cast<std::int64_t>(tensors.size())},
                escalars, tensors, false, kAw3d);
    const AwFamily& f1 = fam[kAw1d];
    if (flags[3])
        section(f1,
                {static_cast<std::int64_t>(f1.Size()), 0,
                 static_cast<std::int64_t>(escalars.size()),
                 static_cast<std::int64_t>(forces.size()), 0},
                escalars, forces, false, kAw1d);
    if (flags[4]) {
        // One subset holding every part, then the material and property tables.
        out.Ints({1});
        out.Text("model", 50);
        out.Ints({0, 0});
        for (const std::size_t slot : {static_cast<std::size_t>(kAw2d), static_cast<std::size_t>(kAw3d),
                                       static_cast<std::size_t>(kAw1d)}) {
            std::vector<std::int64_t> idx;
            for (std::size_t k = 0; k < part_counts[slot]; ++k)
                idx.push_back(static_cast<std::int64_t>(k));
            out.Ints({static_cast<std::int64_t>(idx.size())});
            out.Ints(idx);
        }
        std::set<std::int64_t> mats, props;
        for (const std::size_t slot : {static_cast<std::size_t>(kAw2d), static_cast<std::size_t>(kAw3d),
                                       static_cast<std::size_t>(kAw1d), static_cast<std::size_t>(kAwSph)}) {
            const PartTable pt = part_table(fam[slot]);
            mats.insert(pt.mMats.begin(), pt.mMats.end());
            props.insert(pt.mProps.begin(), pt.mProps.end());
        }
        out.Ints({static_cast<std::int64_t>(mats.size()), static_cast<std::int64_t>(props.size())});
        for (const std::int64_t m : mats)
            out.Text("Material " + std::to_string(m), 50);
        out.Ints(std::vector<std::int64_t>(mats.begin(), mats.end()));
        for (const std::int64_t p : props)
            out.Text("Property " + std::to_string(p), 50);
        out.Ints(std::vector<std::int64_t>(props.begin(), props.end()));
    }
    if (flags[5]) {
        out.Ints({static_cast<std::int64_t>(th_items[0].size()),
                  static_cast<std::int64_t>(th_items[1].size()),
                  static_cast<std::int64_t>(th_items[2].size()),
                  static_cast<std::int64_t>(th_items[3].size())});
        for (const auto& items : th_items) {
            std::vector<std::int64_t> ids;
            for (const auto& item : items)
                ids.push_back(item.first);
            out.Ints(ids);
            for (const auto& item : items)
                out.Text(item.second, 50);
        }
    }
    const AwFamily& fs = fam[kAwSph];
    if (flags[7])
        section(fs,
                {static_cast<std::int64_t>(fs.Size()), 0,
                 static_cast<std::int64_t>(escalars.size()),
                 static_cast<std::int64_t>(tensors.size())},
                escalars, tensors, false, kAwSph);
    return out.Data();
}

void aw_save(const std::string& rPath, const std::string& rBytes) {
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("Could not open file for writing: " + rPath);
    out.write(rBytes.data(), static_cast<std::streamsize>(rBytes.size()));
    out.close();
    if (!out)
        throw WriteError("Radioss animation: failed writing " + rPath);
}

}  // namespace

void write_radioss_anim(const std::string& rPath, const Mesh& rMesh) {
    const detail::ProvenanceSlotlessWrite slotless;
    aw_save(rPath, aw_state_bytes(rMesh, aw_time_of(rMesh)));
}

struct RadiossAnimSeriesWriter::Impl {
    std::filesystem::path mDir;
    std::string mStem;
    std::size_t mWidth = 3;
    std::uint64_t mNext = 1;
    std::size_t mSteps = 0;
};

RadiossAnimSeriesWriter::RadiossAnimSeriesWriter(const std::string& rPath)
    : mpImpl(std::make_unique<Impl>()) {
    const std::filesystem::path path(rPath);
    const std::string name = path.filename().string();
    std::size_t k = name.size();
    while (k > 0 && std::isdigit(static_cast<unsigned char>(name[k - 1])))
        --k;
    mpImpl->mDir = path.parent_path();
    if (name.size() - k < 3 || k < 2 || name[k - 1] != 'A') {
        // Not an animation file name: the name is the stem.
        mpImpl->mStem = name + "A";
        mpImpl->mWidth = 3;
        mpImpl->mNext = 1;
    } else {
        mpImpl->mStem = name.substr(0, k);
        mpImpl->mWidth = name.size() - k;
        mpImpl->mNext = std::stoull(name.substr(k));
    }
}

RadiossAnimSeriesWriter::~RadiossAnimSeriesWriter() = default;
RadiossAnimSeriesWriter::RadiossAnimSeriesWriter(RadiossAnimSeriesWriter&&) noexcept = default;
RadiossAnimSeriesWriter& RadiossAnimSeriesWriter::operator=(RadiossAnimSeriesWriter&&) noexcept =
    default;

void RadiossAnimSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        throw WriteError("Radioss animation: series: the writer was moved from");
    Impl& s = *mpImpl;
    std::string number = std::to_string(s.mNext++);
    if (number.size() < s.mWidth)
        number.insert(0, s.mWidth - number.size(), '0');
    const detail::ProvenanceSlotlessWrite slotless;
    aw_save((s.mDir / (s.mStem + number)).string(), aw_state_bytes(rMesh, Time));
    ++s.mSteps;
}

std::size_t RadiossAnimSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps : 0;
}

void RadiossAnimSeriesWriter::Finalize() {}

}  // namespace meshioplusplus
