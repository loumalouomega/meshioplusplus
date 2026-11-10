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

// System includes
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

// Project includes
#include "tet_remesh.hpp"

// Project includes (private, not installed)
#include "feature_snap.hpp"
#include "predicates.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kTetrmPi = 3.14159265358979323846;

// The outward face of a positive tetrahedron opposite each corner
// (detail/cell_faces' tetra table, indexed by the corner it omits).
constexpr int kTetrmOppositeFace[4][3] = {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 1}};

std::int64_t tetrm_find(std::vector<std::int64_t>& rParent, std::int64_t X) {
    while (rParent[static_cast<std::size_t>(X)] != X) {
        rParent[static_cast<std::size_t>(X)] =
            rParent[static_cast<std::size_t>(rParent[static_cast<std::size_t>(X)])];
        X = rParent[static_cast<std::size_t>(X)];
    }
    return X;
}

void tetrm_union(std::vector<std::int64_t>& rParent, std::int64_t A, std::int64_t B) {
    A = tetrm_find(rParent, A);
    B = tetrm_find(rParent, B);
    if (A == B)
        return;
    // the lower root wins, so the partition's roots do not depend on the order
    if (A < B)
        rParent[static_cast<std::size_t>(B)] = A;
    else
        rParent[static_cast<std::size_t>(A)] = B;
}

Vec3 tetrm_unit_normal(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    return vec3_normalize(vec3_cross(vec3_sub(rB, rA), vec3_sub(rC, rA)));
}

}  // namespace

TriKey tri_key(std::int64_t A, std::int64_t B, std::int64_t C) {
    TriKey k{A, B, C};
    if (k[0] > k[1])
        std::swap(k[0], k[1]);
    if (k[1] > k[2])
        std::swap(k[1], k[2]);
    if (k[0] > k[1])
        std::swap(k[0], k[1]);
    return k;
}

std::size_t TriKeyHash::operator()(const TriKey& rKey) const noexcept {
    std::uint64_t h = 0x9E3779B97F4A7C15ull;
    for (const std::int64_t v : rKey) {
        h ^= static_cast<std::uint64_t>(v) + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
        h *= 0xBF58476D1CE4E5B9ull;
    }
    return static_cast<std::size_t>(h ^ (h >> 31));
}

Vec3 TetRidgeChain::Eval(double S) const {
    if (mPoints.empty())
        return {0.0, 0.0, 0.0};
    if (mPoints.size() == 1 || !(S > 0.0))
        return mPoints.front();
    if (S >= Length())
        return mPoints.back();
    const auto it = std::upper_bound(mArc.begin(), mArc.end(), S);
    const std::size_t i = static_cast<std::size_t>(it - mArc.begin());  // mArc[i-1] <= S < mArc[i]
    const double seg = mArc[i] - mArc[i - 1];
    const double t = seg > 0.0 ? (S - mArc[i - 1]) / seg : 0.0;
    return vec3_add(mPoints[i - 1], vec3_scale(vec3_sub(mPoints[i], mPoints[i - 1]), t));
}

int TetRemeshStore::OrientSign(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    // Shewchuk's orient3d is positive when D lies below the plane of the
    // counter-clockwise (A, B, C): the opposite of meshio's positive volume.
    const double o = orient3d(rA.data(), rB.data(), rC.data(), rD.data());
    return o < 0.0 ? 1 : (o > 0.0 ? -1 : 0);
}

bool TetRemeshStore::Positive(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    return OrientSign(rA, rB, rC, rD) > 0;
}

void TetRemeshStore::BuildVertexTets() {
    const std::size_t n = NumVertices();
    std::vector<std::size_t> count(n, 0);
    for (std::size_t t = 0; t < mTets.size(); ++t)
        if (mTetAlive[t])
            for (const std::int64_t v : mTets[t])
                ++count[static_cast<std::size_t>(v)];
    mVertexTets.assign(n, {});
    for (std::size_t v = 0; v < n; ++v)
        mVertexTets[v].reserve(count[v]);
    for (std::size_t t = 0; t < mTets.size(); ++t)
        if (mTetAlive[t])
            for (const std::int64_t v : mTets[t])
                mVertexTets[static_cast<std::size_t>(v)].push_back(static_cast<std::int64_t>(t));
}

void TetRemeshStore::Build(std::vector<double> Xyz, std::vector<TetCorners> Tets,
                           std::vector<std::int32_t> Class,
                           const std::unordered_map<TriKey, std::int32_t, TriKeyHash>& rFaceTags) {
    const std::size_t n = Xyz.size() / 3;
    if (n >= (std::size_t{1} << 32))
        throw std::invalid_argument("adapt_volume: at most 2^32 points are supported");
    mXyz = std::move(Xyz);
    mTets = std::move(Tets);
    mClass = std::move(Class);
    mTetAlive.assign(mTets.size(), 1);
    for (TetCorners& t : mTets) {
        for (const std::int64_t v : t)
            if (v < 0 || static_cast<std::size_t>(v) >= n)
                throw std::invalid_argument(
                    "adapt_volume: a tetrahedron names a point out of range");
        if (OrientSign(Point(t[0]), Point(t[1]), Point(t[2]), Point(t[3])) < 0)
            std::swap(t[2], t[3]);
    }
    mSize.assign(n, 0.0);
    mKind.assign(n, static_cast<std::uint8_t>(TetVertexKind::Interior));
    mFeature.assign(n, -1);
    mChainParam.assign(n, 0.0);
    mFrozen.assign(n, 0);
    mSource.resize(n);
    std::iota(mSource.begin(), mSource.end(), std::int64_t{0});
    mVertexAlive.assign(n, 0);
    BuildVertexTets();
    for (std::size_t v = 0; v < n; ++v)
        mVertexAlive[v] = !mVertexTets[v].empty();

    // Faces: (key, tet) records sorted, so equal faces are adjacent.
    struct FaceRec {
        TriKey mKey;
        std::int64_t mTet;
    };
    std::vector<FaceRec> recs;
    recs.reserve(mTets.size() * 4);
    for (std::size_t t = 0; t < mTets.size(); ++t)
        for (int k = 0; k < 4; ++k) {
            const auto& f = kTetrmOppositeFace[k];
            recs.push_back({tri_key(mTets[t][f[0]], mTets[t][f[1]], mTets[t][f[2]]),
                            static_cast<std::int64_t>(t)});
        }
    std::sort(recs.begin(), recs.end(), [](const FaceRec& rA, const FaceRec& rB) {
        return rA.mKey != rB.mKey ? rA.mKey < rB.mKey : rA.mTet < rB.mTet;
    });
    mBoundary.clear();
    mBoundary.reserve(recs.size() / 8 + 16);
    for (std::size_t i = 0; i < recs.size();) {
        std::size_t j = i + 1;
        while (j < recs.size() && recs[j].mKey == recs[i].mKey)
            ++j;
        if (j - i > 2)
            throw std::invalid_argument(
                "adapt_volume: a face is shared by more than two tetrahedra (a non-manifold "
                "volume); run check_tetrahedra first");
        TetBoundaryFace face;
        bool boundary = false;
        if (j - i == 1) {
            face.mClassA = mClass[static_cast<std::size_t>(recs[i].mTet)];
            boundary = true;
        } else {
            // two tetrahedra on a face lie on its two sides (exact), unless one
            // is flat; on one side they overlap: a tangled mesh
            const TriKey& k = recs[i].mKey;
            auto apex = [&](std::int64_t T) {
                for (const std::int64_t v : mTets[static_cast<std::size_t>(T)])
                    if (v != k[0] && v != k[1] && v != k[2])
                        return v;
                return k[0];
            };
            const int s1 =
                OrientSign(Point(k[0]), Point(k[1]), Point(k[2]), Point(apex(recs[i].mTet)));
            const int s2 =
                OrientSign(Point(k[0]), Point(k[1]), Point(k[2]), Point(apex(recs[i + 1].mTet)));
            if (s1 * s2 > 0)
                throw std::invalid_argument(
                    "adapt_volume: two tetrahedra overlap across a face (a tangled mesh); run "
                    "untangle first");
            const std::int32_t ca = mClass[static_cast<std::size_t>(recs[i].mTet)];
            const std::int32_t cb = mClass[static_cast<std::size_t>(recs[i + 1].mTet)];
            if (ca != cb) {
                face.mClassA = std::min(ca, cb);
                face.mClassB = std::max(ca, cb);
                boundary = true;
            }
        }
        if (boundary) {
            const auto it = rFaceTags.find(recs[i].mKey);
            if (it != rFaceTags.end())
                face.mTag = it->second;
            mBoundary.emplace(recs[i].mKey, face);
        }
        i = j;
    }
    mRidgeEdges.clear();
    mChains.clear();
    mNumPatches = 0;
}

std::int64_t TetRemeshStore::NumAliveTets() const {
    std::int64_t c = 0;
    for (const std::uint8_t a : mTetAlive)
        c += a;
    return c;
}

void TetRemeshStore::SharedTets(std::int64_t A, std::int64_t B,
                                std::vector<std::int64_t>& rOut) const {
    rOut.clear();
    const auto& la = mVertexTets[static_cast<std::size_t>(A)];
    const auto& lb = mVertexTets[static_cast<std::size_t>(B)];
    std::set_intersection(la.begin(), la.end(), lb.begin(), lb.end(), std::back_inserter(rOut));
}

void TetRemeshStore::VertexRing(std::int64_t V, std::vector<std::int64_t>& rOut) const {
    rOut.clear();
    for (const std::int64_t t : mVertexTets[static_cast<std::size_t>(V)])
        for (const std::int64_t c : mTets[static_cast<std::size_t>(t)])
            if (c != V)
                rOut.push_back(c);
    std::sort(rOut.begin(), rOut.end());
    rOut.erase(std::unique(rOut.begin(), rOut.end()), rOut.end());
}

void TetRemeshStore::BoundaryFacesOf(std::int64_t V, std::vector<TriKey>& rOut) const {
    rOut.clear();
    if (mBoundary.empty())
        return;
    for (const std::int64_t t : mVertexTets[static_cast<std::size_t>(V)]) {
        const TetCorners& c = mTets[static_cast<std::size_t>(t)];
        for (int k = 0; k < 4; ++k) {
            if (c[k] == V)
                continue;  // the face opposite V does not contain it
            const auto& f = kTetrmOppositeFace[k];
            const TriKey key = tri_key(c[f[0]], c[f[1]], c[f[2]]);
            if (mBoundary.count(key))
                rOut.push_back(key);
        }
    }
    std::sort(rOut.begin(), rOut.end());
    rOut.erase(std::unique(rOut.begin(), rOut.end()), rOut.end());
}

bool TetRemeshStore::OnBoundary(std::int64_t V) const {
    if (mBoundary.empty())
        return false;
    for (const std::int64_t t : mVertexTets[static_cast<std::size_t>(V)]) {
        const TetCorners& c = mTets[static_cast<std::size_t>(t)];
        for (int k = 0; k < 4; ++k) {
            if (c[k] == V)
                continue;
            const auto& f = kTetrmOppositeFace[k];
            if (mBoundary.count(tri_key(c[f[0]], c[f[1]], c[f[2]])))
                return true;
        }
    }
    return false;
}

const TetBoundaryFace* TetRemeshStore::FindFace(std::int64_t A, std::int64_t B,
                                                std::int64_t C) const {
    const auto it = mBoundary.find(tri_key(A, B, C));
    return it == mBoundary.end() ? nullptr : &it->second;
}

std::int32_t TetRemeshStore::RidgeChain(std::int64_t A, std::int64_t B) const {
    const auto it = mRidgeEdges.find(edge_key64(A, B));
    return it == mRidgeEdges.end() ? -1 : it->second;
}

double TetRemeshStore::ChainParam(std::int64_t V, std::int32_t Chain, double Near) const {
    const std::size_t v = static_cast<std::size_t>(V);
    if (Kind(V) == TetVertexKind::Ridge && mFeature[v] == Chain)
        return mChainParam[v];
    const TetRidgeChain& c = mChains[static_cast<std::size_t>(Chain)];
    if (V == c.mFront && V != c.mBack)
        return 0.0;
    if (V == c.mBack && V != c.mFront)
        return c.Length();
    return Near < 0.5 * c.Length() ? 0.0 : c.Length();
}

std::array<std::int64_t, 3> TetRemeshStore::OrientedFace(const TriKey& rKey) const {
    const TetBoundaryFace* face = FindFace(rKey[0], rKey[1], rKey[2]);
    std::vector<std::int64_t> ab;
    SharedTets(rKey[0], rKey[1], ab);
    for (const std::int64_t t : ab) {
        const TetCorners& c = mTets[static_cast<std::size_t>(t)];
        int opposite = -1;
        bool has_c = false;
        for (int k = 0; k < 4; ++k) {
            if (c[k] == rKey[2])
                has_c = true;
            else if (c[k] != rKey[0] && c[k] != rKey[1])
                opposite = k;
        }
        if (!has_c || opposite < 0)
            continue;
        if (face && mClass[static_cast<std::size_t>(t)] != face->mClassA)
            continue;
        const auto& f = kTetrmOppositeFace[opposite];
        return {c[f[0]], c[f[1]], c[f[2]]};
    }
    return {rKey[0], rKey[1], rKey[2]};
}

std::int64_t TetRemeshStore::AddVertex(const Vec3& rP, double Size, TetVertexKind Kind,
                                       std::int32_t Feature, double ChainParam,
                                       std::int64_t Source) {
    const std::int64_t id = static_cast<std::int64_t>(NumVertices());
    if (static_cast<std::uint64_t>(id) >= (std::uint64_t{1} << 32))
        throw std::invalid_argument("adapt_volume: at most 2^32 points are supported");
    mXyz.insert(mXyz.end(), {rP[0], rP[1], rP[2]});
    mSize.push_back(Size);
    mKind.push_back(static_cast<std::uint8_t>(Kind));
    mFeature.push_back(Feature);
    mChainParam.push_back(ChainParam);
    mFrozen.push_back(0);
    mVertexAlive.push_back(1);
    mSource.push_back(Source);
    mVertexTets.emplace_back();
    return id;
}

std::int64_t TetRemeshStore::AddTet(const TetCorners& rCorners, std::int32_t Class) {
    const std::int64_t id = static_cast<std::int64_t>(mTets.size());
    mTets.push_back(rCorners);
    mClass.push_back(Class);
    mTetAlive.push_back(1);
    for (const std::int64_t v : rCorners)
        mVertexTets[static_cast<std::size_t>(v)].push_back(id);  // the largest id: stays sorted
    return id;
}

void TetRemeshStore::KillTet(std::int64_t T) {
    const std::size_t t = static_cast<std::size_t>(T);
    if (!mTetAlive[t])
        return;
    mTetAlive[t] = 0;
    for (const std::int64_t v : mTets[t]) {
        auto& l = mVertexTets[static_cast<std::size_t>(v)];
        const auto it = std::lower_bound(l.begin(), l.end(), T);
        if (it != l.end() && *it == T)
            l.erase(it);
    }
}

void TetRemeshStore::ReplaceTets(std::vector<TetCorners> Tets, std::vector<std::int32_t> Class) {
    mTets = std::move(Tets);
    mClass = std::move(Class);
    mTetAlive.assign(mTets.size(), 1);
    BuildVertexTets();
}

TetClassifyReport TetRemeshStore::Classify(double FeatureAngleDeg, double CornerAngleDeg,
                                           bool PreserveFeatures) {
    TetClassifyReport report;
    const std::size_t n = NumVertices();
    std::vector<TriKey> keys;
    keys.reserve(mBoundary.size());
    for (const auto& kv : mBoundary)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    const std::size_t nf = keys.size();

    std::vector<Vec3> normal(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        const auto o = OrientedFace(keys[f]);
        normal[f] = tetrm_unit_normal(Point(o[0]), Point(o[1]), Point(o[2]));
    }

    // Boundary edges and their faces.
    struct EdgeRec {
        std::uint64_t mKey;
        std::int64_t mFace;
    };
    std::vector<EdgeRec> erecs;
    erecs.reserve(nf * 3);
    for (std::size_t f = 0; f < nf; ++f)
        for (int e = 0; e < 3; ++e)
            erecs.push_back(
                {edge_key64(keys[f][e], keys[f][(e + 1) % 3]), static_cast<std::int64_t>(f)});
    std::sort(erecs.begin(), erecs.end(), [](const EdgeRec& rA, const EdgeRec& rB) {
        return rA.mKey != rB.mKey ? rA.mKey < rB.mKey : rA.mFace < rB.mFace;
    });

    const double cos_feature = std::cos(FeatureAngleDeg * kTetrmPi / 180.0);
    std::vector<std::int64_t> parent(nf);
    std::iota(parent.begin(), parent.end(), std::int64_t{0});
    std::vector<std::array<std::int64_t, 2>> ridge_edges;
    std::vector<std::vector<std::int64_t>> ridge_faces;  // the faces of each ridge edge
    for (std::size_t i = 0; i < erecs.size();) {
        std::size_t j = i + 1;
        while (j < erecs.size() && erecs[j].mKey == erecs[i].mKey)
            ++j;
        bool ridge = j - i != 2;
        if (!ridge) {
            const TetBoundaryFace& f1 =
                mBoundary.at(keys[static_cast<std::size_t>(erecs[i].mFace)]);
            const TetBoundaryFace& f2 =
                mBoundary.at(keys[static_cast<std::size_t>(erecs[i + 1].mFace)]);
            ridge = f1.mTag != f2.mTag || f1.mClassA != f2.mClassA || f1.mClassB != f2.mClassB;
            if (!ridge && PreserveFeatures) {
                const Vec3& n1 = normal[static_cast<std::size_t>(erecs[i].mFace)];
                const Vec3& n2 = normal[static_cast<std::size_t>(erecs[i + 1].mFace)];
                const bool valid = vec3_norm_sq(n1) > 0.5 && vec3_norm_sq(n2) > 0.5;
                ridge = valid && vec3_dot(n1, n2) < cos_feature;
            }
        }
        if (ridge) {
            const std::int64_t lo = static_cast<std::int64_t>(erecs[i].mKey >> 32);
            const std::int64_t hi = static_cast<std::int64_t>(erecs[i].mKey & 0xFFFFFFFFull);
            ridge_edges.push_back({lo, hi});
            std::vector<std::int64_t> fs;
            for (std::size_t k = i; k < j; ++k)
                fs.push_back(erecs[k].mFace);
            ridge_faces.push_back(std::move(fs));
        } else {
            tetrm_union(parent, erecs[i].mFace, erecs[i + 1].mFace);
        }
        i = j;
    }

    // Patches, numbered by their first face in key order.
    std::vector<std::int32_t> patch_of_root(nf, -1);
    std::vector<std::int32_t> patch(nf, -1);
    std::int32_t num_patches = 0;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r =
            static_cast<std::size_t>(tetrm_find(parent, static_cast<std::int64_t>(f)));
        if (patch_of_root[r] < 0)
            patch_of_root[r] = num_patches++;
        patch[f] = patch_of_root[r];
        mBoundary.at(keys[f]).mPatch = patch[f];
    }
    mNumPatches = num_patches;

    // Per vertex: the patches of its boundary faces, and its first one.
    std::vector<std::int32_t> first_patch(n, -1);
    std::vector<std::uint8_t> multi_patch(n, 0);
    for (std::size_t f = 0; f < nf; ++f)
        for (const std::int64_t v : keys[f]) {
            std::int32_t& fp = first_patch[static_cast<std::size_t>(v)];
            if (fp < 0)
                fp = patch[f];
            else if (fp != patch[f])
                multi_patch[static_cast<std::size_t>(v)] = 1;
        }

    // The ridge graph: corners and chains, then extra corners where two
    // ridges separating different patches meet at a degree-2 vertex.
    std::vector<Vec3> pts(n);
    for (std::size_t v = 0; v < n; ++v)
        pts[v] = Point(static_cast<std::int64_t>(v));
    const FeatureGraph graph = feature_graph(pts, ridge_edges, CornerAngleDeg);
    std::vector<std::uint8_t> corner(n, 0);
    for (const std::int64_t c : graph.mCorners)
        corner[static_cast<std::size_t>(c)] = 1;
    {
        std::unordered_map<std::uint64_t, std::vector<std::int32_t>> signature;
        signature.reserve(ridge_edges.size() * 2 + 1);
        for (std::size_t e = 0; e < ridge_edges.size(); ++e) {
            std::vector<std::int32_t> sig;
            for (const std::int64_t f : ridge_faces[e])
                sig.push_back(patch[static_cast<std::size_t>(f)]);
            std::sort(sig.begin(), sig.end());
            sig.erase(std::unique(sig.begin(), sig.end()), sig.end());
            signature[edge_key64(ridge_edges[e][0], ridge_edges[e][1])] = std::move(sig);
        }
        for (const auto& chain : graph.mChains)
            for (std::size_t i = 1; i + 1 < chain.size(); ++i)
                if (signature[edge_key64(chain[i - 1], chain[i])] !=
                    signature[edge_key64(chain[i], chain[i + 1])])
                    corner[static_cast<std::size_t>(chain[i])] = 1;
    }
    // A vertex where patches meet without a ridge (a pinch) is a corner too.
    std::vector<std::uint8_t> on_ridge(n, 0);
    for (const auto& e : ridge_edges) {
        on_ridge[static_cast<std::size_t>(e[0])] = 1;
        on_ridge[static_cast<std::size_t>(e[1])] = 1;
    }
    for (std::size_t v = 0; v < n; ++v)
        if (multi_patch[v] && !on_ridge[v])
            corner[v] = 1;

    // Chains, split at every corner.
    mChains.clear();
    mRidgeEdges.clear();
    for (std::size_t v = 0; v < n; ++v) {
        mKind[v] = static_cast<std::uint8_t>(first_patch[v] >= 0 ? TetVertexKind::Surface
                                                                 : TetVertexKind::Interior);
        mFeature[v] = first_patch[v];
        mChainParam[v] = 0.0;
    }
    for (const auto& chain : graph.mChains) {
        std::size_t start = 0;
        for (std::size_t i = 1; i < chain.size(); ++i) {
            if (!corner[static_cast<std::size_t>(chain[i])] && i + 1 < chain.size())
                continue;
            TetRidgeChain rc;
            rc.mFront = chain[start];
            rc.mBack = chain[i];
            corner[static_cast<std::size_t>(rc.mFront)] = 1;  // every chain ends at corners
            corner[static_cast<std::size_t>(rc.mBack)] = 1;
            const std::int32_t id = static_cast<std::int32_t>(mChains.size());
            double arc = 0.0;
            for (std::size_t k = start; k <= i; ++k) {
                const Vec3 p = Point(chain[k]);
                if (k > start)
                    arc += vec3_norm(vec3_sub(p, rc.mPoints.back()));
                rc.mPoints.push_back(p);
                rc.mArc.push_back(arc);
                if (k > start)
                    mRidgeEdges[edge_key64(chain[k - 1], chain[k])] = id;
                if (k > start && k < i) {
                    const std::size_t v = static_cast<std::size_t>(chain[k]);
                    mKind[v] = static_cast<std::uint8_t>(TetVertexKind::Ridge);
                    mFeature[v] = id;
                    mChainParam[v] = arc;
                }
            }
            mChains.push_back(std::move(rc));
            start = i;
        }
    }
    for (std::size_t v = 0; v < n; ++v)
        if (corner[v]) {
            mKind[v] = static_cast<std::uint8_t>(TetVertexKind::Corner);
            mFeature[v] = -1;
        }

    for (std::size_t v = 0; v < n; ++v) {
        if (!mVertexAlive[v])
            continue;
        report.mNumCorners += Kind(static_cast<std::int64_t>(v)) == TetVertexKind::Corner;
        report.mNumRidgeVertices += Kind(static_cast<std::int64_t>(v)) == TetVertexKind::Ridge;
    }
    report.mNumChains = static_cast<std::int64_t>(mChains.size());
    report.mNumPatches = num_patches;
    return report;
}

bool TetRemeshStore::CheckInvariants(std::string& rWhy,
                                     const std::vector<std::uint8_t>& rAllowFlat) const {
    const std::size_t n = NumVertices();
    std::vector<std::size_t> uses(n, 0);
    std::vector<TetCorners> sorted_tets;
    for (std::size_t t = 0; t < mTets.size(); ++t) {
        if (!mTetAlive[t])
            continue;
        const TetCorners& c = mTets[t];
        for (const std::int64_t v : c) {
            if (v < 0 || static_cast<std::size_t>(v) >= n ||
                !mVertexAlive[static_cast<std::size_t>(v)]) {
                rWhy = "tetrahedron " + std::to_string(t) + " names a dead or missing vertex";
                return false;
            }
            const auto& l = mVertexTets[static_cast<std::size_t>(v)];
            if (!std::binary_search(l.begin(), l.end(), static_cast<std::int64_t>(t))) {
                rWhy = "tetrahedron " + std::to_string(t) + " missing from its corner's list";
                return false;
            }
            ++uses[static_cast<std::size_t>(v)];
        }
        const int s = OrientSign(Point(c[0]), Point(c[1]), Point(c[2]), Point(c[3]));
        const bool flat_ok = t < rAllowFlat.size() && rAllowFlat[t];
        if (s < 0 || (s == 0 && !flat_ok)) {
            rWhy = "tetrahedron " + std::to_string(t) + (s < 0 ? " is inverted" : " is flat");
            return false;
        }
        TetCorners k = c;
        std::sort(k.begin(), k.end());
        if (k[0] == k[1] || k[1] == k[2] || k[2] == k[3]) {
            rWhy = "tetrahedron " + std::to_string(t) + " repeats a vertex";
            return false;
        }
        sorted_tets.push_back(k);
    }
    for (std::size_t v = 0; v < n; ++v)
        if (uses[v] != mVertexTets[v].size()) {
            rWhy = "vertex " + std::to_string(v) + " lists a dead tetrahedron";
            return false;
        }
    std::sort(sorted_tets.begin(), sorted_tets.end());
    if (std::adjacent_find(sorted_tets.begin(), sorted_tets.end()) != sorted_tets.end()) {
        rWhy = "two tetrahedra share all four vertices";
        return false;
    }

    struct FaceRec {
        TriKey mKey;
        std::int32_t mClass;
    };
    std::vector<FaceRec> recs;
    for (std::size_t t = 0; t < mTets.size(); ++t)
        if (mTetAlive[t])
            for (int k = 0; k < 4; ++k) {
                const auto& f = kTetrmOppositeFace[k];
                recs.push_back(
                    {tri_key(mTets[t][f[0]], mTets[t][f[1]], mTets[t][f[2]]), mClass[t]});
            }
    std::sort(recs.begin(), recs.end(),
              [](const FaceRec& rA, const FaceRec& rB) { return rA.mKey < rB.mKey; });
    std::size_t num_boundary = 0;
    std::vector<std::uint64_t> outer_edges;
    for (std::size_t i = 0; i < recs.size();) {
        std::size_t j = i + 1;
        while (j < recs.size() && recs[j].mKey == recs[i].mKey)
            ++j;
        if (j - i > 2) {
            rWhy = "a face is shared by more than two tetrahedra";
            return false;
        }
        const bool boundary = j - i == 1 || recs[i].mClass != recs[i + 1].mClass;
        const bool in_map = mBoundary.count(recs[i].mKey) > 0;
        if (boundary != in_map) {
            rWhy = boundary ? "a boundary face is missing from the boundary map"
                            : "the boundary map holds an interior face";
            return false;
        }
        if (boundary)
            ++num_boundary;
        if (j - i == 1) {
            const TriKey& k = recs[i].mKey;
            outer_edges.push_back(edge_key64(k[0], k[1]));
            outer_edges.push_back(edge_key64(k[1], k[2]));
            outer_edges.push_back(edge_key64(k[0], k[2]));
        }
        i = j;
    }
    if (num_boundary != mBoundary.size()) {
        rWhy = "the boundary map holds a face no tetrahedron has";
        return false;
    }
    std::sort(outer_edges.begin(), outer_edges.end());
    for (std::size_t i = 0; i < outer_edges.size();) {
        std::size_t j = i + 1;
        while (j < outer_edges.size() && outer_edges[j] == outer_edges[i])
            ++j;
        if (j - i == 1) {
            rWhy = "the outer boundary has an open edge";
            return false;
        }
        i = j;
    }
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus
