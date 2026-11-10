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
//  Attribution:     Implemented from published descriptions only:
//                   Dapogny, Dobrzynski and Frey, "Three-dimensional adaptive
//                   domain remeshing, implicit domain meshing, and
//                   applications to free and moving boundary problems",
//                   J. Comput. Phys. 262 (2014) (metric edge length, the
//                   split/collapse/swap/move cycle, Hausdorff-bounded surface
//                   changes); Klingner and Shewchuk, "Aggressive Tetrahedral
//                   Mesh Improvement", IMR 2007, and Klingner's thesis
//                   UCB/EECS-2008-145 (schedule and quality measures);
//                   Freitag and Ollivier-Gooch, IJNME 40(21), 1997; Hu et al.,
//                   "Fast Tetrahedral Meshing in the Wild", ACM TOG 39(4),
//                   2020 (conformal AMIPS); Jiang et al., "Declarative
//                   Specification for Unstructured Mesh Editing Algorithms",
//                   ACM TOG 41(6), 2022 (precondition / edit / postcondition
//                   / rollback); Botsch and Kobbelt 2004 (the 4/5 and 4/3
//                   length thresholds); Alliez et al., "Variational
//                   Tetrahedral Meshing", SIGGRAPH 2005 (ODT). No code of a
//                   copyleft or non-commercial mesher (MMG, CGAL, TetGen,
//                   Gmsh, Netgen, Mesquite, TetWild, Stellar) is read or
//                   vendored.

// System includes
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/adapt_volume.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "adapt_volume_internal.hpp"
#include "../detail/decimate_shared.hpp"
#include "../detail/fem_kernels.hpp"
#include "../detail/size_field_utils.hpp"
#include "../detail/tet_flips.hpp"
#include "../detail/tet_remesh.hpp"
#include "../detail/tet_scope.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

using detail::SizeFieldUtils;
using detail::TetCorners;
using detail::TetRemeshStore;
using detail::TetVertexKind;
using detail::TriKey;
using detail::Vec3;

constexpr const char* kAdvPrefix = "meshio++: adapt_volume: ";

/// The integer cell_data names read as materials when the caller names none.
constexpr const char* kAdvDefaultMaterials[] = {"octree_mesh:material", "medit:ref",
                                                "gmsh:physical", "gmsh:geometrical", "cell_tags"};

/// The local facet (cell_faces' tetra order) of a tetrahedron opposite each corner.
constexpr int kAdvFacetOpposite[4] = {1, 2, 0, 3};
/// The corners of each local facet of a tetrahedron (cell_faces' tetra table).
constexpr int kAdvFacetCorners[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};

using AdvClock = std::chrono::steady_clock;

double adv_seconds(AdvClock::time_point Since) {
    return std::chrono::duration<double>(AdvClock::now() - Since).count();
}

/// A material class: the input block, the cell regions and the material values
/// a tetrahedron carries, all inherited by every tetrahedron made from it.
struct AdvClass {
    std::size_t mBlock = 0;
    std::vector<std::int32_t> mRegions;    ///< cell-region indices, ascending
    std::vector<std::int64_t> mMaterials;  ///< one value per material array
    std::int64_t mRepresentative = 0;      ///< a local input cell of mBlock with these values
    std::string mLabel;
};

/// The measures the report takes before and after.
struct AdvMeasures {
    double mMinDihedral = 0.0;
    double mMaxDihedral = 0.0;
    double mMeanQuality = 0.0;
    double mVolume = 0.0;
    std::vector<double> mClassVolume;
};

/// Why a candidate edit was refused (the first guard that refused it).
enum class AdvReject : std::uint8_t {
    None = 0,
    NotEdge,
    Link,
    Inversion,
    Quality,
    Length,
    Class,
    Hausdorff,
    Normal,
};

/// A candidate edge in a pass's queue: stale once either end's version moved.
struct AdvEdgeCand {
    double mKey = 0.0;  ///< the metric length (collapse) or its negative (split)
    std::int64_t mA = 0;
    std::int64_t mB = 0;
    std::uint32_t mVersionA = 0;
    std::uint32_t mVersionB = 0;

    bool operator>(const AdvEdgeCand& rO) const {
        if (mKey != rO.mKey)
            return mKey > rO.mKey;
        if (mA != rO.mA)
            return mA > rO.mA;
        return mB > rO.mB;
    }
};

using AdvQueue =
    std::priority_queue<AdvEdgeCand, std::vector<AdvEdgeCand>, std::greater<AdvEdgeCand>>;

/// Edges shorter than this (in the metric) are collapsed, longer than kAdvLong split.
constexpr double kAdvShort = 0.8;
constexpr double kAdvLong = 4.0 / 3.0;

/**
 * The remesher: the store, the size field, the classes and tags, the data
 * carried per vertex, and one method per pass.
 */
class AdvRemesher {
public:
    AdvRemesher(const Mesh& rMesh, const AdaptVolumeOptions& rOptions, AdaptVolumeResult& rResult)
        : mrMesh(rMesh), mrOptions(rOptions), mrResult(rResult) {}

    /// Setup, classification, sizes and the per-vertex state; idempotent.
    void Prepare();
    void Run();
    /// The guard refusing the collapse A -> B (test hook).
    AdvReject ProbeCollapse(std::int64_t A, std::int64_t B) { return TryCollapse(A, B); }

private:
    const Mesh& mrMesh;
    const AdaptVolumeOptions& mrOptions;
    AdaptVolumeResult& mrResult;

    TetRemeshStore mStore;
    std::vector<AdvClass> mClasses;
    std::vector<std::string> mMaterialNames;
    /// Per tag id: the (side region, owner class) pairs a boundary face carries.
    std::vector<std::vector<std::pair<std::int32_t, std::int32_t>>> mTags;
    /// Per point region: membership of every store vertex.
    std::vector<std::vector<std::uint8_t>> mPointMember;
    std::vector<std::size_t> mPointRegions;
    std::vector<detail::DecimFloatData> mFloatData;
    std::size_t mNumInputPoints = 0;
    double mDiagonal = 1.0;
    double mHausd = 0.0;
    double mQualityMin = 0.0;
    double mEpsVol6 = 0.0;
    std::map<std::string, double> mTimes;
    std::vector<std::uint32_t> mVersion;  ///< per vertex, bumped by every edit of its star
    bool mCheck = false;  ///< MESHIOPLUSPLUS_ADAPT_VOLUME_CHECK: verify after each pass
    std::vector<TetCorners> mFlatKeys;  ///< input tetrahedra that were already flat, sorted

    // scratch, reused across candidates
    std::vector<std::int64_t> mScratchShared, mScratchRingA, mScratchRingB, mScratchTetsA;
    std::vector<std::int64_t> mScratchCommon, mScratchLinkV;
    std::vector<std::uint64_t> mScratchEdgesA, mScratchEdgesB, mScratchLinkE, mScratchEdgeCommon;
    std::vector<TriKey> mScratchTriA, mScratchTriB, mScratchFaces;

    void Setup();
    void ReadClasses(std::vector<std::int32_t>& rClassOfCell);
    void ReadTags(const std::vector<std::int32_t>& rClassOfCell,
                  const std::vector<TetCorners>& rTets,
                  std::unordered_map<TriKey, std::int32_t, detail::TriKeyHash>& rFaceTags);
    void BuildSizes();
    AdvMeasures Measure() const;
    void Emit();
    void Verify(const char* pStage) const;
    void CountReject(AdvReject Why);

    double EdgeLength(std::int64_t A, std::int64_t B) const {
        return SizeFieldUtils::EdgeLength(mStore.Point(A), mStore.Point(B),
                                          mStore.mSize[static_cast<std::size_t>(A)],
                                          mStore.mSize[static_cast<std::size_t>(B)]);
    }
    /// Every edge of an alive tetrahedron, sorted, as (lo, hi).
    std::vector<std::array<std::int64_t, 2>> AliveEdges() const;
    void PushEdge(AdvQueue& rQueue, std::int64_t A, std::int64_t B, double Key) const;

    // collapse
    AdvReject LinkCondition(std::int64_t A, std::int64_t B);
    AdvReject TryCollapse(std::int64_t A, std::int64_t B);
    void CommitCollapse(std::int64_t A, std::int64_t B);
    std::int64_t CollapsePass();

    // split
    /// A new vertex on the edge (a, b) at fraction @p T, carrying interpolated data.
    std::int64_t NewVertex(std::int64_t A, std::int64_t B, double T, const Vec3& rP,
                           TetVertexKind Kind, std::int32_t Feature, double ChainParam);
    /// The boundary faces through the edge (a, b) (the shared tets in mScratchShared).
    void EdgeBoundaryFaces(std::int64_t A, std::int64_t B, std::vector<TriKey>& rOut) const;
    AdvReject TrySplit(std::int64_t A, std::int64_t B, Vec3& rP, double& rT);
    void CommitSplit(std::int64_t A, std::int64_t B, const Vec3& rP, double T);
    std::int64_t SplitPass();

    // swaps and the shared smoothing (detail/tet_flips.hpp, per class)
    /// The alive tetrahedra of each class that a move may touch, and the rest.
    void TetsByClass(std::vector<std::vector<detail::TetFlipTet>>& rByClass,
                     std::vector<detail::TetFlipTet>& rFixed,
                     std::vector<std::int32_t>& rFixedClass) const;
    void StoreTetsByClass(std::vector<std::vector<detail::TetFlipTet>>& rByClass,
                          std::vector<detail::TetFlipTet>& rFixed,
                          std::vector<std::int32_t>& rFixedClass);
    std::int64_t SwapPass();

    void Time(const char* pStage, AdvClock::time_point Since) {
        mTimes[pStage] += adv_seconds(Since);
    }

public:
    /// The guarded quality of a tetrahedron (higher is better); -1 when not positive.
    double Quality(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) const {
        if (detail::tet_volume6(rA, rB, rC, rD) <= mEpsVol6 ||
            !TetRemeshStore::Positive(rA, rB, rC, rD))
            return -1.0;
        if (mrOptions.mQualityMetric == AdaptVolumeQuality::MinDihedral)
            return detail::tet_min_dihedral(rA, rB, rC, rD);
        return detail::tet_mean_ratio(rA, rB, rC, rD);
    }
    double TetQuality(const TetCorners& rT) const {
        return Quality(mStore.Point(rT[0]), mStore.Point(rT[1]), mStore.Point(rT[2]),
                       mStore.Point(rT[3]));
    }
};

// --- setup -----------------------------------------------------------------

void AdvRemesher::ReadClasses(std::vector<std::int32_t>& rClassOfCell) {
    const std::vector<std::int64_t> bases = detail::block_bases(mrMesh);
    const std::size_t total = static_cast<std::size_t>(detail::total_cells(bases));

    // Material arrays: the caller's, or the known names present.
    if (!mrOptions.mMaterialArrays.empty()) {
        for (const std::string& name : mrOptions.mMaterialArrays)
            if (!mrMesh.HasCellData(name))
                throw std::invalid_argument(std::string(kAdvPrefix) + "material array '" + name +
                                            "' is not a cell_data array");
        mMaterialNames = mrOptions.mMaterialArrays;
    } else {
        for (const char* name : kAdvDefaultMaterials)
            if (mrMesh.HasCellData(name))
                mMaterialNames.emplace_back(name);
    }
    for (const std::string& name : mMaterialNames) {
        if (mrMesh.CellDataNumBlocks(name) != mrMesh.NumCellBlocks())
            throw std::invalid_argument(std::string(kAdvPrefix) + "material array '" + name +
                                        "' does not have one entry per block");
        for (std::size_t b = 0; b < mrMesh.NumCellBlocks(); ++b) {
            const NDArray& a = mrMesh.CellData(name, b);
            if (detail::is_float_dtype(a.Dtype()) ||
                detail::rows(a) != mrMesh.Cells(b).NumCells() || (a.Size() != detail::rows(a)))
                throw std::invalid_argument(std::string(kAdvPrefix) + "material array '" + name +
                                            "' must be one integer per cell");
        }
    }

    // Cell-region memberships.
    std::vector<std::vector<std::int32_t>> regions_of(total);
    for (std::size_t r = 0; r < mrMesh.NumRegions(); ++r) {
        const Region& reg = mrMesh.Region(r);
        if (reg.mKind != RegionKind::Cell)
            continue;
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i)
            if (e[i] >= 0 && static_cast<std::size_t>(e[i]) < total)
                regions_of[static_cast<std::size_t>(e[i])].push_back(static_cast<std::int32_t>(r));
    }

    rClassOfCell.assign(total, 0);
    std::map<std::vector<std::int64_t>, std::int32_t> ids;
    for (std::size_t b = 0; b < mrMesh.NumCellBlocks(); ++b) {
        const std::size_t nc = mrMesh.Cells(b).NumCells();
        std::vector<const NDArray*> mats;
        for (const std::string& name : mMaterialNames)
            mats.push_back(&mrMesh.CellData(name, b));
        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t g = static_cast<std::size_t>(bases[b]) + c;
            std::vector<std::int64_t> key{static_cast<std::int64_t>(b)};
            for (const std::int32_t r : regions_of[g])
                key.push_back(r);
            key.push_back(-1);
            std::vector<std::int64_t> values;
            for (const NDArray* a : mats)
                values.push_back(detail::read_int(*a, c));
            key.insert(key.end(), values.begin(), values.end());
            const auto it = ids.find(key);
            if (it != ids.end()) {
                rClassOfCell[g] = it->second;
                continue;
            }
            const std::int32_t id = static_cast<std::int32_t>(mClasses.size());
            ids.emplace(key, id);
            AdvClass cls;
            cls.mBlock = b;
            cls.mRegions = regions_of[g];
            cls.mMaterials = values;
            cls.mRepresentative = static_cast<std::int64_t>(c);
            cls.mLabel = "block " + std::to_string(b);
            for (const std::int32_t r : cls.mRegions)
                cls.mLabel += " " + mrMesh.Region(static_cast<std::size_t>(r)).mName;
            for (std::size_t m = 0; m < values.size(); ++m)
                cls.mLabel += " " + mMaterialNames[m] + "=" + std::to_string(values[m]);
            mClasses.push_back(std::move(cls));
            rClassOfCell[g] = id;
        }
    }
}

void AdvRemesher::ReadTags(
    const std::vector<std::int32_t>& rClassOfCell, const std::vector<TetCorners>& rTets,
    std::unordered_map<TriKey, std::int32_t, detail::TriKeyHash>& rFaceTags) {
    // Per face: the (side region, owner class) pairs, then one id per distinct list.
    std::map<TriKey, std::vector<std::pair<std::int32_t, std::int32_t>>> lists;
    for (std::size_t r = 0; r < mrMesh.NumRegions(); ++r) {
        const Region& reg = mrMesh.Region(r);
        if (reg.mKind != RegionKind::Side)
            continue;
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i) {
            const std::int64_t cell = e[2 * i], facet = e[2 * i + 1];
            if (cell < 0 || static_cast<std::size_t>(cell) >= rTets.size() || facet < 0 ||
                facet > 3)
                continue;
            const TetCorners& t = rTets[static_cast<std::size_t>(cell)];
            const auto& f = kAdvFacetCorners[facet];
            lists[detail::tri_key(t[f[0]], t[f[1]], t[f[2]])].push_back(
                {static_cast<std::int32_t>(r), rClassOfCell[static_cast<std::size_t>(cell)]});
        }
    }
    std::map<std::vector<std::pair<std::int32_t, std::int32_t>>, std::int32_t> ids;
    mTags.assign(1, {});  // tag 0: none
    for (auto& [key, list] : lists) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
        auto it = ids.find(list);
        if (it == ids.end()) {
            it = ids.emplace(list, static_cast<std::int32_t>(mTags.size())).first;
            mTags.push_back(list);
        }
        rFaceTags.emplace(key, it->second);
    }
}

void AdvRemesher::Setup() {
    detail::require_tet_only_blocks(mrMesh, kAdvPrefix);
    if (mrMesh.PointDim() != 3)
        throw std::invalid_argument(std::string(kAdvPrefix) + "needs 3-D points");
    if (!(mrOptions.mHgrad >= 1.0))
        throw std::invalid_argument(std::string(kAdvPrefix) + "hgrad must be >= 1");
    if (mrOptions.mMaxIterations < 0)
        throw std::invalid_argument(std::string(kAdvPrefix) + "max_iterations must be >= 0");
    if (!(mrOptions.mMaxNormalDeviation > 0.0) || mrOptions.mMaxNormalDeviation > 180.0)
        throw std::invalid_argument(std::string(kAdvPrefix) +
                                    "max_normal_deviation must be in (0, 180] degrees");
    const std::size_t n = mrMesh.NumPoints();
    mNumInputPoints = n;
    if (!mrOptions.mFrozen.empty() && mrOptions.mFrozen.size() != n)
        throw std::invalid_argument(std::string(kAdvPrefix) +
                                    "frozen must have one entry per point");

    std::vector<double> xyz(n * 3);
    {
        const detail::DoubleView pv(mrMesh.Points());
        for (std::size_t i = 0; i < n * 3; ++i)
            xyz[i] = pv[i];
    }
    std::vector<TetCorners> tets;
    {
        const std::vector<std::int64_t> bases = detail::block_bases(mrMesh);
        tets.resize(static_cast<std::size_t>(detail::total_cells(bases)));
        std::size_t bi = 0;
        for (const auto cb : mrMesh.CellRange()) {
            const detail::Int64View cv(cb.Conn());
            const std::size_t base = static_cast<std::size_t>(bases[bi++]);
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                for (std::size_t k = 0; k < 4; ++k)
                    tets[base + c][k] = cv[c * 4 + k];
        }
    }
    std::vector<std::int32_t> class_of_cell;
    ReadClasses(class_of_cell);
    std::unordered_map<TriKey, std::int32_t, detail::TriKeyHash> face_tags;
    ReadTags(class_of_cell, tets, face_tags);

    // Scale.
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    for (std::size_t i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d) {
            const double x = xyz[i * 3 + static_cast<std::size_t>(d)];
            lo[d] = i == 0 ? x : std::min(lo[d], x);
            hi[d] = i == 0 ? x : std::max(hi[d], x);
        }
    mDiagonal = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                          (hi[2] - lo[2]) * (hi[2] - lo[2]));
    if (!(mDiagonal > 0.0))
        mDiagonal = 1.0;
    mEpsVol6 = mDiagonal * mDiagonal * mDiagonal * 1e-14;
    if (mrOptions.mHausd > 0.0)
        mHausd = mrOptions.mHausdRelative ? mrOptions.mHausd * mDiagonal : mrOptions.mHausd;
    else
        mHausd = 1e-3 * mDiagonal;
    if (mrOptions.mQualityMin >= 0.0)
        mQualityMin = mrOptions.mQualityMin;
    else
        mQualityMin = mrOptions.mQualityMetric == AdaptVolumeQuality::MinDihedral ? 10.0 : 0.2;

    mStore.Build(std::move(xyz), std::move(tets), std::move(class_of_cell), face_tags);
    const detail::TetClassifyReport cr = mStore.Classify(
        mrOptions.mFeatureAngle, mrOptions.mCornerAngle, mrOptions.mPreserveFeatures);
    mrResult.mNumCorners = cr.mNumCorners;
    mrResult.mNumRidgeVertices = cr.mNumRidgeVertices;
    mrResult.mNumChains = cr.mNumChains;
    mrResult.mNumPatches = cr.mNumPatches;

    for (std::size_t v = 0; v < n; ++v) {
        const bool frozen = !mrOptions.mFrozen.empty() && mrOptions.mFrozen[v];
        const bool boundary =
            mrOptions.mPreserveBoundary && mStore.OnBoundary(static_cast<std::int64_t>(v));
        mStore.mFrozen[v] = static_cast<std::uint8_t>(frozen || boundary);
    }

    // Data carried per vertex.
    mFloatData = detail::decim_float_working_copies(mrMesh, n);
    for (std::size_t r = 0; r < mrMesh.NumRegions(); ++r) {
        const Region& reg = mrMesh.Region(r);
        if (reg.mKind != RegionKind::Point)
            continue;
        std::vector<std::uint8_t> member(n, 0);
        const std::int64_t* e = reg.Entries();
        for (std::size_t i = 0; i < reg.NumEntries(); ++i)
            if (e[i] >= 0 && static_cast<std::size_t>(e[i]) < n)
                member[static_cast<std::size_t>(e[i])] = 1;
        mPointRegions.push_back(r);
        mPointMember.push_back(std::move(member));
    }
    BuildSizes();
}

void AdvRemesher::BuildSizes() {
    const std::size_t n = mStore.NumVertices();
    std::vector<double> boundary_sum(n, 0.0);
    std::vector<std::int32_t> boundary_count(n, 0);
    std::vector<TriKey> keys;
    keys.reserve(mStore.mBoundary.size());
    for (const auto& kv : mStore.mBoundary)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    std::vector<std::uint64_t> bedges;
    for (const TriKey& k : keys)
        for (int e = 0; e < 3; ++e)
            bedges.push_back(detail::edge_key64(k[e], k[(e + 1) % 3]));
    std::sort(bedges.begin(), bedges.end());
    bedges.erase(std::unique(bedges.begin(), bedges.end()), bedges.end());
    for (const std::uint64_t e : bedges) {
        const std::int64_t a = static_cast<std::int64_t>(e >> 32);
        const std::int64_t b = static_cast<std::int64_t>(e & 0xFFFFFFFFull);
        const double l = detail::vec3_norm(detail::vec3_sub(mStore.Point(a), mStore.Point(b)));
        boundary_sum[static_cast<std::size_t>(a)] += l;
        boundary_sum[static_cast<std::size_t>(b)] += l;
        ++boundary_count[static_cast<std::size_t>(a)];
        ++boundary_count[static_cast<std::size_t>(b)];
    }

    // The caller's sizes.
    std::vector<double> user;
    if (!mrOptions.mSizeArray.empty()) {
        if (!mrMesh.HasPointData(mrOptions.mSizeArray))
            throw std::invalid_argument(std::string(kAdvPrefix) + "size_array '" +
                                        mrOptions.mSizeArray + "' is not a point_data array");
        const NDArray& a = mrMesh.PointData(mrOptions.mSizeArray);
        if (a.Size() != mNumInputPoints)
            throw std::invalid_argument(std::string(kAdvPrefix) + "size_array '" +
                                        mrOptions.mSizeArray + "' must be one value per point");
        user.resize(mNumInputPoints);
        for (std::size_t i = 0; i < mNumInputPoints; ++i)
            user[i] = detail::read_double(a, i);
    }

    // Curvature: the turn of the vertex normals along each boundary edge of a patch.
    std::vector<double> curvature(n, 0.0);
    if (mrOptions.mCurvatureSizing) {
        std::vector<Vec3> vnormal(n, Vec3{0.0, 0.0, 0.0});
        for (const TriKey& k : keys) {
            const auto o = mStore.OrientedFace(k);
            const Vec3 a = mStore.Point(o[0]), b = mStore.Point(o[1]), c = mStore.Point(o[2]);
            const Vec3 an = detail::vec3_cross(detail::vec3_sub(b, a), detail::vec3_sub(c, a));
            const std::int32_t patch = mStore.mBoundary.at(k).mPatch;
            for (const std::int64_t v : k)
                if (mStore.Kind(v) == TetVertexKind::Surface &&
                    mStore.mFeature[static_cast<std::size_t>(v)] == patch)
                    vnormal[static_cast<std::size_t>(v)] =
                        detail::vec3_add(vnormal[static_cast<std::size_t>(v)], an);
        }
        for (Vec3& v : vnormal)
            v = detail::vec3_normalize(v);
        for (const std::uint64_t e : bedges) {
            const std::int64_t a = static_cast<std::int64_t>(e >> 32);
            const std::int64_t b = static_cast<std::int64_t>(e & 0xFFFFFFFFull);
            const Vec3& na = vnormal[static_cast<std::size_t>(a)];
            const Vec3& nb = vnormal[static_cast<std::size_t>(b)];
            if (detail::vec3_norm_sq(na) < 0.5 || detail::vec3_norm_sq(nb) < 0.5)
                continue;
            const double l = detail::vec3_norm(detail::vec3_sub(mStore.Point(a), mStore.Point(b)));
            if (!(l > 0.0))
                continue;
            const double k = detail::vec3_norm(detail::vec3_sub(na, nb)) / l;
            curvature[static_cast<std::size_t>(a)] =
                std::max(curvature[static_cast<std::size_t>(a)], k);
            curvature[static_cast<std::size_t>(b)] =
                std::max(curvature[static_cast<std::size_t>(b)], k);
        }
    }

    // Seeds.
    double hmin_auto = std::numeric_limits<double>::infinity();
    std::vector<double> seed(n, 0.0);
    for (std::size_t v = 0; v < n; ++v)
        if (boundary_count[v] > 0) {
            double h = boundary_sum[v] / boundary_count[v];
            if (!user.empty() && v < user.size() && user[v] > 0.0)
                h = std::min(h, user[v]);
            if (curvature[v] > 0.0)
                h = std::min(h, std::sqrt(8.0 * mHausd / curvature[v]));
            seed[v] = h;
            hmin_auto = std::min(hmin_auto, h);
        }
    const double hmax = mrOptions.mHmax > 0.0 ? mrOptions.mHmax : mDiagonal / 8.0;
    double hmin = mrOptions.mHmin > 0.0 ? mrOptions.mHmin
                                        : (std::isfinite(hmin_auto) ? hmin_auto : hmax / 32.0);
    hmin = std::min(hmin, hmax);
    if (mrOptions.mHmin > 0.0 && mrOptions.mHmax > 0.0 && mrOptions.mHmin > mrOptions.mHmax)
        throw std::invalid_argument(std::string(kAdvPrefix) + "hmin must not exceed hmax");
    for (std::size_t v = 0; v < n; ++v)
        if (boundary_count[v] == 0) {
            double h = hmax;
            if (!user.empty() && v < user.size() && user[v] > 0.0)
                h = std::min(h, user[v]);
            seed[v] = h;
        }

    // Grade over the tetrahedra's edges, then clamp.
    std::vector<std::array<std::int64_t, 2>> edges;
    edges.reserve(mStore.NumTets() * 6);
    for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
        const TetCorners& c = mStore.mTets[t];
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                edges.push_back({std::min(c[i], c[j]), std::max(c[i], c[j])});
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    SizeFieldUtils::Grade(mStore.mXyz, edges, mrOptions.mHgrad, seed);
    for (double& h : seed)
        h = std::clamp(h, hmin, hmax);
    mStore.mSize = std::move(seed);
    mrResult.mHmin = hmin;
    mrResult.mHmax = hmax;
    mrResult.mHausd = mHausd;
}

AdvMeasures AdvRemesher::Measure() const {
    AdvMeasures m;
    m.mClassVolume.assign(mClasses.size(), 0.0);
    std::vector<std::size_t> alive;
    for (std::size_t t = 0; t < mStore.NumTets(); ++t)
        if (mStore.mTetAlive[t])
            alive.push_back(t);
    std::vector<double> mn(alive.size()), mx(alive.size()), q(alive.size()), vol(alive.size());
    parallel_for(alive.size(), [&](std::size_t i) {
        const TetCorners& c = mStore.mTets[alive[i]];
        const Vec3 a = mStore.Point(c[0]), b = mStore.Point(c[1]), cc = mStore.Point(c[2]),
                   d = mStore.Point(c[3]);
        detail::tet_dihedral_range(a, b, cc, d, mn[i], mx[i]);
        q[i] = mrOptions.mQualityMetric == AdaptVolumeQuality::MinDihedral
                   ? detail::tet_min_dihedral(a, b, cc, d)
                   : detail::tet_mean_ratio(a, b, cc, d);
        vol[i] = detail::tet_volume(a, b, cc, d);
    });
    m.mMinDihedral = alive.empty() ? 0.0 : 180.0;
    double qsum = 0.0;
    for (std::size_t i = 0; i < alive.size(); ++i) {
        m.mMinDihedral = std::min(m.mMinDihedral, mn[i]);
        m.mMaxDihedral = std::max(m.mMaxDihedral, mx[i]);
        qsum += q[i];
        m.mVolume += vol[i];
        m.mClassVolume[static_cast<std::size_t>(mStore.mClass[alive[i]])] += vol[i];
    }
    m.mMeanQuality = alive.empty() ? 0.0 : qsum / static_cast<double>(alive.size());
    return m;
}

// --- shared pass machinery ----------------------------------------------------

void AdvRemesher::Verify(const char* pStage) const {
    if (!mCheck)
        return;
    std::vector<std::uint8_t> allow_flat(mStore.NumTets(), 0);
    if (!mFlatKeys.empty())
        for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
            TetCorners k = mStore.mTets[t];
            std::sort(k.begin(), k.end());
            allow_flat[t] = std::binary_search(mFlatKeys.begin(), mFlatKeys.end(), k);
        }
    std::string why;
    if (!mStore.CheckInvariants(why, allow_flat))
        throw std::logic_error(std::string(kAdvPrefix) + "invariant broken after " + pStage + ": " +
                               why);
}

void AdvRemesher::CountReject(AdvReject Why) {
    switch (Why) {
        case AdvReject::Link:
            ++mrResult.mNumRejectedByLink;
            break;
        case AdvReject::Inversion:
            ++mrResult.mNumRejectedByInversion;
            break;
        case AdvReject::Quality:
            ++mrResult.mNumRejectedByQuality;
            break;
        case AdvReject::Length:
            ++mrResult.mNumRejectedByLength;
            break;
        case AdvReject::Class:
            ++mrResult.mNumRejectedByClass;
            break;
        case AdvReject::Hausdorff:
            ++mrResult.mNumRejectedByHausdorff;
            break;
        case AdvReject::Normal:
            ++mrResult.mNumRejectedByNormal;
            break;
        default:
            break;
    }
}

std::vector<std::array<std::int64_t, 2>> AdvRemesher::AliveEdges() const {
    std::vector<std::array<std::int64_t, 2>> edges;
    edges.reserve(static_cast<std::size_t>(mStore.NumAliveTets()) * 6);
    for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
        if (!mStore.mTetAlive[t])
            continue;
        const TetCorners& c = mStore.mTets[t];
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                edges.push_back({std::min(c[i], c[j]), std::max(c[i], c[j])});
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

void AdvRemesher::PushEdge(AdvQueue& rQueue, std::int64_t A, std::int64_t B, double Key) const {
    const std::int64_t lo = std::min(A, B), hi = std::max(A, B);
    rQueue.push({Key, lo, hi, mVersion[static_cast<std::size_t>(lo)],
                 mVersion[static_cast<std::size_t>(hi)]});
}

// --- collapse ------------------------------------------------------------------

// The link condition of the edge (a, b) (Dey et al.; Edelsbrunner): the
// vertices, edges and triangles common to the links of a and b are exactly
// those of the link of ab. A triangle common to both links would leave two
// tetrahedra on the same four vertices (the duplicate-tet check).
AdvReject AdvRemesher::LinkCondition(std::int64_t A, std::int64_t B) {
    const auto& shared = mScratchShared;
    // vertices
    mStore.VertexRing(A, mScratchRingA);
    mStore.VertexRing(B, mScratchRingB);
    mScratchCommon.clear();
    std::set_intersection(mScratchRingA.begin(), mScratchRingA.end(), mScratchRingB.begin(),
                          mScratchRingB.end(), std::back_inserter(mScratchCommon));
    mScratchLinkV.clear();
    mScratchLinkE.clear();
    for (const std::int64_t t : shared) {
        std::int64_t o[2];
        int k = 0;
        for (const std::int64_t c : mStore.mTets[static_cast<std::size_t>(t)])
            if (c != A && c != B)
                o[k++] = c;
        mScratchLinkV.push_back(o[0]);
        mScratchLinkV.push_back(o[1]);
        mScratchLinkE.push_back(detail::edge_key64(o[0], o[1]));
    }
    std::sort(mScratchLinkV.begin(), mScratchLinkV.end());
    mScratchLinkV.erase(std::unique(mScratchLinkV.begin(), mScratchLinkV.end()),
                        mScratchLinkV.end());
    if (mScratchCommon != mScratchLinkV)
        return AdvReject::Link;
    // edges and triangles of the links
    auto link_of = [&](std::int64_t V, std::int64_t Other, std::vector<std::uint64_t>& rEdges,
                       std::vector<TriKey>& rTris) {
        rEdges.clear();
        rTris.clear();
        for (const std::int64_t t : mStore.mVertexTets[static_cast<std::size_t>(V)]) {
            std::int64_t o[3];
            int k = 0;
            bool has_other = false;
            for (const std::int64_t c : mStore.mTets[static_cast<std::size_t>(t)]) {
                if (c == V)
                    continue;
                has_other = has_other || c == Other;
                o[k++] = c;
            }
            rEdges.push_back(detail::edge_key64(o[0], o[1]));
            rEdges.push_back(detail::edge_key64(o[1], o[2]));
            rEdges.push_back(detail::edge_key64(o[0], o[2]));
            if (!has_other)
                rTris.push_back(detail::tri_key(o[0], o[1], o[2]));
        }
        std::sort(rEdges.begin(), rEdges.end());
        rEdges.erase(std::unique(rEdges.begin(), rEdges.end()), rEdges.end());
        std::sort(rTris.begin(), rTris.end());
    };
    link_of(A, B, mScratchEdgesA, mScratchTriA);
    link_of(B, A, mScratchEdgesB, mScratchTriB);
    mScratchEdgeCommon.clear();
    std::set_intersection(mScratchEdgesA.begin(), mScratchEdgesA.end(), mScratchEdgesB.begin(),
                          mScratchEdgesB.end(), std::back_inserter(mScratchEdgeCommon));
    std::sort(mScratchLinkE.begin(), mScratchLinkE.end());
    for (const std::uint64_t e : mScratchEdgeCommon)
        if (!std::binary_search(mScratchLinkE.begin(), mScratchLinkE.end(), e))
            return AdvReject::Link;
    // a triangle opposite a and opposite b: the collapse would duplicate a tetrahedron
    std::size_t i = 0, j = 0;
    while (i < mScratchTriA.size() && j < mScratchTriB.size()) {
        if (mScratchTriA[i] == mScratchTriB[j])
            return AdvReject::Link;
        if (mScratchTriA[i] < mScratchTriB[j])
            ++i;
        else
            ++j;
    }
    return AdvReject::None;
}

// Precondition, simulated edit and postconditions of the half-edge collapse
// a -> b (a removed, its tetrahedra moved onto b); nothing is changed.
AdvReject AdvRemesher::TryCollapse(std::int64_t A, std::int64_t B) {
    const std::size_t a = static_cast<std::size_t>(A);
    if (mStore.mFrozen[a] || mStore.Kind(A) != TetVertexKind::Interior)
        return AdvReject::Class;
    mStore.SharedTets(A, B, mScratchShared);
    if (mScratchShared.empty())
        return AdvReject::NotEdge;
    if (const AdvReject r = LinkCondition(A, B); r != AdvReject::None)
        return r;
    // the tetrahedra that move, and their orientation (exact)
    mScratchTetsA.clear();
    double old_min = std::numeric_limits<double>::infinity();
    for (const std::int64_t t : mStore.mVertexTets[a]) {
        old_min = std::min(old_min, TetQuality(mStore.mTets[static_cast<std::size_t>(t)]));
        if (!std::binary_search(mScratchShared.begin(), mScratchShared.end(), t))
            mScratchTetsA.push_back(t);
    }
    double new_min = std::numeric_limits<double>::infinity();
    for (const std::int64_t t : mScratchTetsA) {
        TetCorners c = mStore.mTets[static_cast<std::size_t>(t)];
        for (std::int64_t& v : c)
            if (v == A)
                v = B;
        const Vec3 p0 = mStore.Point(c[0]), p1 = mStore.Point(c[1]), p2 = mStore.Point(c[2]),
                   p3 = mStore.Point(c[3]);
        if (detail::tet_volume6(p0, p1, p2, p3) <= mEpsVol6 ||
            !TetRemeshStore::Positive(p0, p1, p2, p3))
            return AdvReject::Inversion;
        new_min = std::min(new_min, Quality(p0, p1, p2, p3));
    }
    // no new edge longer than the split threshold
    for (const std::int64_t x : mScratchRingA) {
        if (x == B || std::binary_search(mScratchRingB.begin(), mScratchRingB.end(), x))
            continue;
        if (EdgeLength(B, x) > kAdvLong)
            return AdvReject::Length;
    }
    if (new_min < std::min(old_min, mQualityMin))
        return AdvReject::Quality;
    return AdvReject::None;
}

void AdvRemesher::CommitCollapse(std::int64_t A, std::int64_t B) {
    const std::size_t a = static_cast<std::size_t>(A);
    mStore.SharedTets(A, B, mScratchShared);
    for (const std::int64_t t : mScratchShared)
        mStore.KillTet(t);
    const std::vector<std::int64_t> moving = mStore.mVertexTets[a];
    auto& lb = mStore.mVertexTets[static_cast<std::size_t>(B)];
    for (const std::int64_t t : moving) {
        for (std::int64_t& v : mStore.mTets[static_cast<std::size_t>(t)])
            if (v == A)
                v = B;
        detail::decim_sorted_insert(lb, t);
    }
    mStore.mVertexTets[a].clear();
    mStore.mVertexAlive[a] = 0;
    ++mVersion[a];
    ++mVersion[static_cast<std::size_t>(B)];
}

std::int64_t AdvRemesher::CollapsePass() {
    AdvQueue queue;
    {
        const auto edges = AliveEdges();
        std::vector<double> len(edges.size());
        parallel_for(edges.size(),
                     [&](std::size_t i) { len[i] = EdgeLength(edges[i][0], edges[i][1]); });
        for (std::size_t i = 0; i < edges.size(); ++i)
            if (len[i] < kAdvShort)
                PushEdge(queue, edges[i][0], edges[i][1], len[i]);
    }
    std::int64_t applied = 0;
    std::vector<std::int64_t> ring;
    while (!queue.empty()) {
        const AdvEdgeCand e = queue.top();
        queue.pop();
        const std::size_t ea = static_cast<std::size_t>(e.mA), eb = static_cast<std::size_t>(e.mB);
        if (!mStore.mVertexAlive[ea] || !mStore.mVertexAlive[eb] || mVersion[ea] != e.mVersionA ||
            mVersion[eb] != e.mVersionB)
            continue;  // stale
        // remove the less constrained end first (the lower kind survives), then the other
        std::int64_t first = e.mB, second = e.mA;
        if (mStore.mKind[ea] > mStore.mKind[eb])
            std::swap(first, second);
        AdvReject why = TryCollapse(first, second);
        std::int64_t removed = first, kept = second;
        if (why != AdvReject::None) {
            const AdvReject other = TryCollapse(second, first);
            if (other == AdvReject::None) {
                why = other;
                removed = second;
                kept = first;
            }
        }
        if (why == AdvReject::NotEdge)
            continue;
        if (why != AdvReject::None) {
            CountReject(why);
            continue;
        }
        CommitCollapse(removed, kept);
        ++applied;
        mStore.VertexRing(kept, ring);
        for (const std::int64_t x : ring) {
            const double l = EdgeLength(kept, x);
            if (l < kAdvShort)
                PushEdge(queue, kept, x, l);
        }
    }
    return applied;
}

// --- split ---------------------------------------------------------------------

std::int64_t AdvRemesher::NewVertex(std::int64_t A, std::int64_t B, double T, const Vec3& rP,
                                    TetVertexKind Kind, std::int32_t Feature, double ChainParam) {
    const std::size_t a = static_cast<std::size_t>(A), b = static_cast<std::size_t>(B);
    const double h = SizeFieldUtils::Interpolate(mStore.mSize[a], mStore.mSize[b], T);
    const std::int64_t source = T <= 0.5 ? mStore.mSource[a] : mStore.mSource[b];
    const std::int64_t m = mStore.AddVertex(rP, h, Kind, Feature, ChainParam, source);
    for (detail::DecimFloatData& f : mFloatData) {
        const std::size_t nc = f.mNumComponents;
        for (std::size_t k = 0; k < nc; ++k)
            f.mValues.push_back((1.0 - T) * f.mValues[a * nc + k] + T * f.mValues[b * nc + k]);
    }
    // A point region takes the new vertex when it holds both ends and the new
    // vertex has the kind of one of them (so the midpoint of an edge between
    // two corners does not join a region of corners).
    const std::uint8_t ka = mStore.mKind[a], kb = mStore.mKind[b];
    const std::uint8_t km = static_cast<std::uint8_t>(Kind);
    for (auto& member : mPointMember)
        member.push_back(
            static_cast<std::uint8_t>(member[a] && member[b] && (km == ka || km == kb)));
    mVersion.push_back(0);
    return m;
}

void AdvRemesher::EdgeBoundaryFaces(std::int64_t A, std::int64_t B,
                                    std::vector<TriKey>& rOut) const {
    rOut.clear();
    if (mStore.mBoundary.empty())
        return;
    for (const std::int64_t t : mScratchShared)
        for (const std::int64_t x : mStore.mTets[static_cast<std::size_t>(t)]) {
            if (x == A || x == B)
                continue;
            const TriKey k = detail::tri_key(A, B, x);
            if (mStore.mBoundary.count(k))
                rOut.push_back(k);
        }
    std::sort(rOut.begin(), rOut.end());
    rOut.erase(std::unique(rOut.begin(), rOut.end()), rOut.end());
}

// Precondition, simulated edit and postconditions of splitting the edge (a, b)
// at its metric midpoint: every child positive and no worse than the floor.
AdvReject AdvRemesher::TrySplit(std::int64_t A, std::int64_t B, Vec3& rP, double& rT) {
    mStore.SharedTets(A, B, mScratchShared);
    if (mScratchShared.empty())
        return AdvReject::NotEdge;
    EdgeBoundaryFaces(A, B, mScratchFaces);
    if (!mScratchFaces.empty() && mrOptions.mPreserveBoundary)
        return AdvReject::Class;
    const std::size_t a = static_cast<std::size_t>(A), b = static_cast<std::size_t>(B);
    rT = SizeFieldUtils::MidpointParameter(mStore.mSize[a], mStore.mSize[b]);
    const Vec3 pa = mStore.Point(A), pb = mStore.Point(B);
    rP = detail::vec3_add(pa, detail::vec3_scale(detail::vec3_sub(pb, pa), rT));
    double old_min = std::numeric_limits<double>::infinity();
    double new_min = std::numeric_limits<double>::infinity();
    for (const std::int64_t t : mScratchShared) {
        const TetCorners& c = mStore.mTets[static_cast<std::size_t>(t)];
        old_min = std::min(old_min, TetQuality(c));
        for (const std::int64_t replaced : {B, A}) {
            Vec3 q[4];
            for (int k = 0; k < 4; ++k)
                q[k] = c[k] == replaced ? rP : mStore.Point(c[k]);
            if (detail::tet_volume6(q[0], q[1], q[2], q[3]) <= mEpsVol6 ||
                !TetRemeshStore::Positive(q[0], q[1], q[2], q[3]))
                return AdvReject::Inversion;
            new_min = std::min(new_min, Quality(q[0], q[1], q[2], q[3]));
        }
    }
    if (new_min < std::min(old_min, mQualityMin))
        return AdvReject::Quality;
    return AdvReject::None;
}

void AdvRemesher::CommitSplit(std::int64_t A, std::int64_t B, const Vec3& rP, double T) {
    mStore.SharedTets(A, B, mScratchShared);
    EdgeBoundaryFaces(A, B, mScratchFaces);
    const std::size_t a = static_cast<std::size_t>(A), b = static_cast<std::size_t>(B);
    TetVertexKind kind = TetVertexKind::Interior;
    std::int32_t feature = -1;
    double param = 0.0;
    const std::int32_t chain = mStore.RidgeChain(A, B);
    if (chain >= 0) {
        kind = TetVertexKind::Ridge;
        feature = chain;
        const double sa = mStore.ChainParam(A, chain, mStore.ChainParam(B, chain, 0.0));
        const double sb = mStore.ChainParam(B, chain, sa);
        param = sa + (sb - sa) * T;
    } else if (!mScratchFaces.empty()) {
        kind = TetVertexKind::Surface;
        feature = mStore.mBoundary.at(mScratchFaces.front()).mPatch;
    }
    const std::int64_t m = NewVertex(A, B, T, rP, kind, feature, param);
    const std::vector<std::int64_t> shared = mScratchShared;
    for (const std::int64_t t : shared) {
        const TetCorners c = mStore.mTets[static_cast<std::size_t>(t)];
        const std::int32_t cls = mStore.mClass[static_cast<std::size_t>(t)];
        TetCorners lo = c, hi = c;
        for (int k = 0; k < 4; ++k) {
            if (c[k] == B)
                lo[k] = m;
            if (c[k] == A)
                hi[k] = m;
        }
        mStore.KillTet(t);
        mStore.AddTet(lo, cls);
        mStore.AddTet(hi, cls);
    }
    for (const TriKey& k : mScratchFaces) {
        const detail::TetBoundaryFace face = mStore.mBoundary.at(k);
        std::int64_t x = k[0];
        for (const std::int64_t v : k)
            if (v != A && v != B)
                x = v;
        mStore.mBoundary.erase(k);
        mStore.mBoundary.emplace(detail::tri_key(A, m, x), face);
        mStore.mBoundary.emplace(detail::tri_key(m, B, x), face);
    }
    if (chain >= 0) {
        mStore.mRidgeEdges.erase(detail::edge_key64(A, B));
        mStore.mRidgeEdges[detail::edge_key64(A, m)] = chain;
        mStore.mRidgeEdges[detail::edge_key64(m, B)] = chain;
    }
    ++mVersion[a];
    ++mVersion[b];
}

std::int64_t AdvRemesher::SplitPass() {
    AdvQueue queue;  // keyed by the negative length: longest first
    {
        const auto edges = AliveEdges();
        std::vector<double> len(edges.size());
        parallel_for(edges.size(),
                     [&](std::size_t i) { len[i] = EdgeLength(edges[i][0], edges[i][1]); });
        for (std::size_t i = 0; i < edges.size(); ++i)
            if (len[i] > kAdvLong)
                PushEdge(queue, edges[i][0], edges[i][1], -len[i]);
    }
    // a bound on the work of one pass, against a size field far below the mesh
    const std::int64_t budget = 4 * mStore.NumAliveTets() + 1024;
    std::int64_t applied = 0;
    std::vector<std::int64_t> ring;
    while (!queue.empty() && applied < budget) {
        const AdvEdgeCand e = queue.top();
        queue.pop();
        const std::size_t ea = static_cast<std::size_t>(e.mA), eb = static_cast<std::size_t>(e.mB);
        if (!mStore.mVertexAlive[ea] || !mStore.mVertexAlive[eb] || mVersion[ea] != e.mVersionA ||
            mVersion[eb] != e.mVersionB)
            continue;
        Vec3 p;
        double t = 0.5;
        const AdvReject why = TrySplit(e.mA, e.mB, p, t);
        if (why == AdvReject::NotEdge)
            continue;
        if (why != AdvReject::None) {
            CountReject(why);
            continue;
        }
        CommitSplit(e.mA, e.mB, p, t);
        ++applied;
        const std::int64_t m = static_cast<std::int64_t>(mStore.NumVertices()) - 1;
        mStore.VertexRing(m, ring);
        for (const std::int64_t x : ring) {
            const double l = EdgeLength(m, x);
            if (l > kAdvLong)
                PushEdge(queue, m, x, -l);
        }
    }
    return applied;
}

// --- swaps ---------------------------------------------------------------------

/// The least dihedral angle below which a tetrahedron is a target of the
/// moves (optimize_volume's worst-cell threshold).
constexpr double kAdvWorstThresholdDeg = 40.0;

void AdvRemesher::TetsByClass(std::vector<std::vector<detail::TetFlipTet>>& rByClass,
                              std::vector<detail::TetFlipTet>& rFixed,
                              std::vector<std::int32_t>& rFixedClass) const {
    rByClass.assign(mClasses.size(), {});
    rFixed.clear();
    rFixedClass.clear();
    // tagged faces between two tetrahedra of one class: keep both sides out
    bool any_inner = false;
    for (const auto& kv : mStore.mBoundary)
        any_inner = any_inner || kv.second.mClassA == kv.second.mClassB;
    for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
        if (!mStore.mTetAlive[t])
            continue;
        const TetCorners& c = mStore.mTets[t];
        bool fixed = false;
        if (any_inner)
            for (int k = 0; k < 4 && !fixed; ++k) {
                const detail::TetBoundaryFace* f =
                    mStore.FindFace(c[(k + 1) % 4], c[(k + 2) % 4], c[(k + 3) % 4]);
                fixed = f && f->mClassA == f->mClassB;
            }
        if (fixed) {
            rFixed.push_back(c);
            rFixedClass.push_back(mStore.mClass[t]);
        } else {
            rByClass[static_cast<std::size_t>(mStore.mClass[t])].push_back(c);
        }
    }
}

void AdvRemesher::StoreTetsByClass(std::vector<std::vector<detail::TetFlipTet>>& rByClass,
                                   std::vector<detail::TetFlipTet>& rFixed,
                                   std::vector<std::int32_t>& rFixedClass) {
    std::vector<TetCorners> tets;
    std::vector<std::int32_t> cls;
    for (std::size_t c = 0; c < rByClass.size(); ++c) {
        tets.insert(tets.end(), rByClass[c].begin(), rByClass[c].end());
        cls.insert(cls.end(), rByClass[c].size(), static_cast<std::int32_t>(c));
    }
    tets.insert(tets.end(), rFixed.begin(), rFixed.end());
    cls.insert(cls.end(), rFixedClass.begin(), rFixedClass.end());
    mStore.ReplaceTets(std::move(tets), std::move(cls));
    for (std::uint32_t& v : mVersion)
        ++v;
}

// optimize_volume's worst-cell moves, run on each class alone: a face between
// two classes has one tetrahedron in each vector, so it is never flipped, and
// the boundary faces are untouched. Every move raises the least dihedral
// angle of the tetrahedra it replaces.
std::int64_t AdvRemesher::SwapPass() {
    std::vector<std::vector<detail::TetFlipTet>> by_class;
    std::vector<detail::TetFlipTet> fixed;
    std::vector<std::int32_t> fixed_class;
    TetsByClass(by_class, fixed, fixed_class);
    constexpr double kMinImprovement = 1e-6;
    constexpr int kMaxRing = 7;
    const auto quality = [&](const detail::TetFlipTet& t) {
        return detail::tet_flip_worst_q(mStore.mXyz, t, mEpsVol6);
    };
    std::int64_t applied = 0;
    for (auto& tets : by_class) {
        if (tets.empty())
            continue;
        applied += detail::tet_flip_pass_23(mStore.mXyz, tets, kMinImprovement, mEpsVol6, quality,
                                            kAdvWorstThresholdDeg);
        std::int64_t n32 = 0;
        applied += detail::tet_flip_pass_edge_removal(mStore.mXyz, tets, kMinImprovement, mEpsVol6,
                                                      kMaxRing, kAdvWorstThresholdDeg, n32);
        applied += detail::tet_flip_pass_multi_face(mStore.mXyz, tets, kMinImprovement, mEpsVol6,
                                                    kMaxRing, kAdvWorstThresholdDeg);
    }
    if (applied > 0)
        StoreTetsByClass(by_class, fixed, fixed_class);
    return applied;
}

// --- output ----------------------------------------------------------------

void AdvRemesher::Emit() {
    const std::size_t nv = mStore.NumVertices();
    // Vertices: the used ones, in store order (survivors of the input, then new).
    std::vector<std::uint8_t> used(nv, 0);
    for (std::size_t t = 0; t < mStore.NumTets(); ++t)
        if (mStore.mTetAlive[t])
            for (const std::int64_t v : mStore.mTets[t])
                used[static_cast<std::size_t>(v)] = 1;
    std::vector<std::int64_t> remap(nv, -1);
    std::vector<std::int64_t> order;
    for (std::size_t v = 0; v < nv; ++v)
        if (used[v]) {
            remap[v] = static_cast<std::int64_t>(order.size());
            order.push_back(static_cast<std::int64_t>(v));
        }
    const std::size_t np = order.size();
    Mesh& out = mrResult.mMesh;
    {
        const NDArray& points = mrMesh.Points();
        NDArray pts = NDArray::Uninit(points.Dtype(), {np, 3});
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t d = 0; d < 3; ++d)
                detail::write_double(pts, i * 3 + d,
                                     mStore.mXyz[static_cast<std::size_t>(order[i]) * 3 + d]);
        out.AssignPoints(std::move(pts));
    }

    // Tetrahedra: one block per input block, in store order.
    const std::size_t nb = mrMesh.NumCellBlocks();
    std::vector<std::vector<std::int64_t>> block_tets(nb);
    for (std::size_t t = 0; t < mStore.NumTets(); ++t)
        if (mStore.mTetAlive[t])
            block_tets[mClasses[static_cast<std::size_t>(mStore.mClass[t])].mBlock].push_back(
                static_cast<std::int64_t>(t));
    std::vector<std::int64_t> out_index(mStore.NumTets(), -1);
    std::vector<std::int64_t> base(nb + 1, 0);
    for (std::size_t b = 0; b < nb; ++b) {
        base[b + 1] = base[b] + static_cast<std::int64_t>(block_tets[b].size());
        NDArray conn = NDArray::Uninit(DType::Int64, {block_tets[b].size(), 4});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < block_tets[b].size(); ++i) {
            const std::size_t t = static_cast<std::size_t>(block_tets[b][i]);
            out_index[t] = base[b] + static_cast<std::int64_t>(i);
            for (std::size_t k = 0; k < 4; ++k)
                dst[i * 4 + k] = remap[static_cast<std::size_t>(mStore.mTets[t][k])];
        }
        out.AddCellBlock(cell_type_name(CellType::Tetra), std::move(conn));
    }

    // point_data: float arrays carried (interpolated on splits), the others
    // from the input point each vertex keeps.
    std::vector<std::int64_t> source(np);
    for (std::size_t i = 0; i < np; ++i)
        source[i] = mStore.mSource[static_cast<std::size_t>(order[i])];
    for (const std::string& name : mrMesh.PointDataNames()) {
        const NDArray& a = mrMesh.PointData(name);
        const auto fit =
            std::find_if(mFloatData.begin(), mFloatData.end(),
                         [&](const detail::DecimFloatData& rF) { return rF.mName == name; });
        if (fit == mFloatData.end()) {
            if (detail::rows(a) == mNumInputPoints)
                out.AddPointData(name, detail::subset_gather_rows(a, source));
            continue;
        }
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = np;
        NDArray o = NDArray::Uninit(a.Dtype(), shape);
        const std::size_t nc = fit->mNumComponents;
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < nc; ++k)
                detail::write_double(o, i * nc + k,
                                     fit->mValues[static_cast<std::size_t>(order[i]) * nc + k]);
        out.AddPointData(name, std::move(o));
    }

    // cell_data: the material arrays by class; the rest has no correspondence.
    for (const std::string& name : mMaterialNames) {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < nb; ++b) {
            std::vector<std::int64_t> rows;
            rows.reserve(block_tets[b].size());
            for (const std::int64_t t : block_tets[b])
                rows.push_back(
                    mClasses[static_cast<std::size_t>(mStore.mClass[static_cast<std::size_t>(t)])]
                        .mRepresentative);
            blocks.push_back(detail::subset_gather_rows(mrMesh.CellData(name, b), rows));
        }
        out.AddCellData(name, std::move(blocks));
    }
    std::vector<std::string> dropped;
    for (const std::string& name : mrMesh.CellDataNames())
        if (std::find(mMaterialNames.begin(), mMaterialNames.end(), name) == mMaterialNames.end())
            dropped.push_back(name);
    if (!dropped.empty()) {
        std::string msg = "dropped cell_data with no correspondence across the edits:";
        for (const std::string& d : dropped)
            msg += " " + d;
        log::warn("{}{}", kAdvPrefix, msg);
        mrResult.mWarnings.push_back(msg);
    }
    for (const std::string& name : mrMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(mrMesh.FieldData(name)));

    // Regions.
    std::vector<std::vector<std::int64_t>> cell_entries(mrMesh.NumRegions());
    for (std::size_t t = 0; t < mStore.NumTets(); ++t)
        if (mStore.mTetAlive[t])
            for (const std::int32_t r :
                 mClasses[static_cast<std::size_t>(mStore.mClass[t])].mRegions)
                cell_entries[static_cast<std::size_t>(r)].push_back(out_index[t]);
    std::vector<std::vector<std::int64_t>> side_entries(mrMesh.NumRegions());
    {
        std::vector<TriKey> keys;
        for (const auto& kv : mStore.mBoundary)
            if (kv.second.mTag > 0)
                keys.push_back(kv.first);
        std::sort(keys.begin(), keys.end());
        std::vector<std::int64_t> ab;
        for (const TriKey& k : keys) {
            const auto& list = mTags[static_cast<std::size_t>(mStore.mBoundary.at(k).mTag)];
            mStore.SharedTets(k[0], k[1], ab);
            for (const auto& [region, owner] : list) {
                for (const std::int64_t t : ab) {
                    const TetCorners& c = mStore.mTets[static_cast<std::size_t>(t)];
                    int opposite = -1;
                    bool has = false;
                    for (int q = 0; q < 4; ++q) {
                        if (c[q] == k[2])
                            has = true;
                        else if (c[q] != k[0] && c[q] != k[1])
                            opposite = q;
                    }
                    if (!has || opposite < 0 || mStore.mClass[static_cast<std::size_t>(t)] != owner)
                        continue;
                    side_entries[static_cast<std::size_t>(region)].push_back(
                        out_index[static_cast<std::size_t>(t)]);
                    side_entries[static_cast<std::size_t>(region)].push_back(
                        kAdvFacetOpposite[opposite]);
                    break;
                }
            }
        }
    }
    std::size_t point_region = 0;
    for (std::size_t r = 0; r < mrMesh.NumRegions(); ++r) {
        const Region& reg = mrMesh.Region(r);
        std::vector<std::int64_t> entries;
        if (reg.mKind == RegionKind::Cell) {
            entries = std::move(cell_entries[r]);
        } else if (reg.mKind == RegionKind::Side) {
            entries = std::move(side_entries[r]);
        } else {
            const auto& member = mPointMember[point_region++];
            for (std::size_t v = 0; v < np; ++v)
                if (member[static_cast<std::size_t>(order[v])])
                    entries.push_back(static_cast<std::int64_t>(v));
        }
        const std::size_t stride = reg.mKind == RegionKind::Side ? 2 : 1;
        NDArray e = reg.mKind == RegionKind::Side
                        ? NDArray::Uninit(DType::Int64, {entries.size() / stride, 2})
                        : NDArray::Uninit(DType::Int64, {entries.size()});
        if (!entries.empty())
            std::copy(entries.begin(), entries.end(), e.As<std::int64_t>());
        out.AddRegion(Region(reg.mName, reg.mKind, reg.mDim, reg.mTag, std::move(e)));
    }
    for (std::size_t i = 0; i < mrMesh.NumPropertySets(); ++i)
        out.AddPropertySet(mrMesh.GetPropertySet(i));

    mrResult.mNumPointsAfter = static_cast<std::int64_t>(np);
    mrResult.mNumTetsAfter = base[nb];
}

// --- driver ----------------------------------------------------------------

void AdvRemesher::Prepare() {
    if (const char* p_env = std::getenv("MESHIOPLUSPLUS_ADAPT_VOLUME_CHECK"))
        mCheck = p_env[0] != '\0' && p_env[0] != '0';
    Setup();
    mVersion.assign(mStore.NumVertices(), 0);
    mFlatKeys.clear();
    for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
        const TetCorners& c = mStore.mTets[t];
        if (TetRemeshStore::OrientSign(mStore.Point(c[0]), mStore.Point(c[1]), mStore.Point(c[2]),
                                       mStore.Point(c[3])) == 0) {
            TetCorners k = c;
            std::sort(k.begin(), k.end());
            mFlatKeys.push_back(k);
        }
    }
    std::sort(mFlatKeys.begin(), mFlatKeys.end());
    Verify("setup");
}

void AdvRemesher::Run() {
    auto t0 = AdvClock::now();
    Prepare();
    mrResult.mNumTetsBefore = mStore.NumAliveTets();
    {
        std::int64_t used = 0;
        for (const std::uint8_t a : mStore.mVertexAlive)
            used += a;
        mrResult.mNumPointsBefore = used;
    }
    const AdvMeasures before = Measure();
    mrResult.mMinDihedralBefore = before.mMinDihedral;
    mrResult.mMaxDihedralBefore = before.mMaxDihedral;
    mrResult.mMeanQualityBefore = before.mMeanQuality;
    mrResult.mVolumeBefore = before.mVolume;
    Time("setup", t0);

    for (int it = 0; it < mrOptions.mMaxIterations; ++it) {
        std::int64_t changed = 0;
        if (mrOptions.mEnableCollapse) {
            t0 = AdvClock::now();
            const std::int64_t n = CollapsePass();
            mrResult.mNumCollapsed += n;
            changed += n;
            Time("collapse", t0);
            Verify("collapse");
        }
        if (mrOptions.mEnableSplit) {
            t0 = AdvClock::now();
            const std::int64_t n = SplitPass();
            mrResult.mNumSplit += n;
            changed += n;
            Time("split", t0);
            Verify("split");
        }
        if (mrOptions.mEnableSwap) {
            t0 = AdvClock::now();
            const std::int64_t n = SwapPass();
            mrResult.mNumSwapped += n;
            changed += n;
            Time("swap", t0);
            Verify("swap");
        }
        ++mrResult.mIterations;
        if (static_cast<double>(changed) < 0.01 * static_cast<double>(mStore.NumAliveTets()))
            break;
        if (mrOptions.mTargetMinDihedral > 0.0 &&
            Measure().mMinDihedral >= mrOptions.mTargetMinDihedral)
            break;
    }

    t0 = AdvClock::now();
    const AdvMeasures after = Measure();
    mrResult.mMinDihedralAfter = after.mMinDihedral;
    mrResult.mMaxDihedralAfter = after.mMaxDihedral;
    mrResult.mMeanQualityAfter = after.mMeanQuality;
    mrResult.mVolumeAfter = after.mVolume;
    for (std::size_t c = 0; c < mClasses.size(); ++c)
        mrResult.mRegionVolumeChange.emplace_back(mClasses[c].mLabel,
                                                  after.mClassVolume[c] - before.mClassVolume[c]);
    Emit();
    Time("output", t0);
    for (const char* stage : {"setup", "collapse", "split", "swap", "smooth", "polish", "output"})
        mrResult.mTimings.emplace_back(stage, mTimes.count(stage) ? mTimes[stage] : 0.0);
}

}  // namespace

AdaptVolumeQuality adapt_volume_quality_from_name(const std::string& rName) {
    if (rName == "min_dihedral")
        return AdaptVolumeQuality::MinDihedral;
    if (rName == "mean_ratio")
        return AdaptVolumeQuality::MeanRatio;
    if (rName == "amips")
        return AdaptVolumeQuality::Amips;
    throw std::invalid_argument(std::string(kAdvPrefix) +
                                "quality_metric must be 'min_dihedral', "
                                "'mean_ratio' or 'amips', got '" +
                                rName + "'");
}

const char* adapt_volume_quality_name(AdaptVolumeQuality Quality) {
    switch (Quality) {
        case AdaptVolumeQuality::MeanRatio:
            return "mean_ratio";
        case AdaptVolumeQuality::Amips:
            return "amips";
        case AdaptVolumeQuality::MinDihedral:
        default:
            return "min_dihedral";
    }
}

namespace detail {

std::string adapt_volume_probe_collapse(const Mesh& rMesh, const AdaptVolumeOptions& rOptions,
                                        std::int64_t A, std::int64_t B) {
    AdaptVolumeResult result;
    AdvRemesher remesher(rMesh, rOptions, result);
    remesher.Prepare();
    if (A < 0 || B < 0 || A == B || static_cast<std::size_t>(std::max(A, B)) >= rMesh.NumPoints())
        throw std::invalid_argument("adapt_volume_probe_collapse: bad vertex ids");
    switch (remesher.ProbeCollapse(A, B)) {
        case AdvReject::None:
            return "none";
        case AdvReject::NotEdge:
            return "not_edge";
        case AdvReject::Link:
            return "link";
        case AdvReject::Inversion:
            return "inversion";
        case AdvReject::Quality:
            return "quality";
        case AdvReject::Length:
            return "length";
        case AdvReject::Class:
            return "class";
        case AdvReject::Hausdorff:
            return "hausdorff";
        case AdvReject::Normal:
            return "normal";
    }
    return "none";
}

}  // namespace detail

AdaptVolumeResult adapt_volume(const Mesh& rMesh, const AdaptVolumeOptions& rOptions) {
    AdaptVolumeResult result;
    AdvRemesher(rMesh, rOptions, result).Run();
    return result;
}

}  // namespace meshioplusplus
