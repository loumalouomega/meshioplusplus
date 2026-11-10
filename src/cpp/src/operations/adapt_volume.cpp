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
#include <cstdio>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
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
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/smooth.hpp"
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
    // the reference skin: the input's boundary and interfaces, by patch
    detail::TriangleSoup mRefSoup;
    std::unique_ptr<detail::DistanceQuery> mRefQuery;
    std::vector<std::int32_t> mRefPatch;
    std::vector<std::array<std::int64_t, 3>> mRefTriangles;
    double mCosNormal = 0.0;
    /// Reference vertices, edge midpoints and centroids, and a grid over them
    /// (the reverse half of the Hausdorff bound).
    std::vector<Vec3> mRefSamples;
    std::unique_ptr<detail::SpatialGrid> mRefGrid;
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

    // the boundary: reference skin, envelope, normals
    void BuildReference();
    /// The nearest point of patch @p Patch of the reference skin to @p rP.
    bool ProjectToPatch(const Vec3& rP, std::int32_t Patch, Vec3& rOut) const;
    /// Every sample of every triangle within hausd of its patch of the reference.
    bool InEnvelope(const std::vector<std::array<Vec3, 3>>& rTris,
                    const std::vector<std::int32_t>& rPatches) const;
    /// The boundary faces of V (each wound outward from its class-A side) and
    /// the same faces with V at @p rP, their patches, and whether every face's
    /// normal turns by at most the bound.
    /**
     * @brief The reverse half of the bound: every reference sample within
     * hausd of the faces @p rOld stays within max(0.9 hausd, its distance to
     * them) of @p rNew or of the boundary faces around the vertices of
     * @p rAround other than @p rReplaced (sorted), the faces the edit replaces.
     */
    bool ReverseCovered(const std::vector<std::array<Vec3, 3>>& rOld,
                        const std::vector<std::array<Vec3, 3>>& rNew,
                        const std::vector<std::int64_t>& rAround,
                        const std::vector<TriKey>& rReplaced);
    /// The current geometry of the boundary faces @p rKeys.
    std::vector<std::array<Vec3, 3>> FaceGeometry(const std::vector<TriKey>& rKeys) const;
    bool MovedFaces(std::int64_t V, std::int64_t Onto, const Vec3& rP,
                    std::vector<std::array<Vec3, 3>>& rNew, std::vector<std::int32_t>& rPatches);
    AdvReject BoundaryCollapse(std::int64_t A, std::int64_t B);
    void CommitBoundaryCollapse(std::int64_t A, std::int64_t B);
    std::int64_t BoundarySmoothPass();
    /// The star of V stays positive with V at rP and its least quality is at least @p Floor.
    bool StarValid(std::int64_t V, const Vec3& rP, double Floor) const;
    double StarMinQuality(std::int64_t V) const;
    double BoundaryDeviation() const;

    // smoothing
    /// Per class, 0 for the interior vertices of that class a smoothing may move.
    std::vector<std::uint8_t> InteriorPins(std::int32_t Class) const;
    std::int64_t SmoothInteriorPass(bool Odt);
    double StarAmips(std::int64_t V, const Vec3& rP, double& rMinQuality, bool& rValid) const;
    std::int64_t AmipsPass();

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
    const TetVertexKind ka = mStore.Kind(A);
    if (mStore.mFrozen[a] || ka == TetVertexKind::Corner)
        return AdvReject::Class;
    mStore.SharedTets(A, B, mScratchShared);
    if (mScratchShared.empty())
        return AdvReject::NotEdge;
    // a surface vertex moves along an edge of its patch, a ridge vertex along its chain
    if (ka == TetVertexKind::Surface) {
        EdgeBoundaryFaces(A, B, mScratchFaces);
        if (mScratchFaces.size() != 2)
            return AdvReject::Class;
        for (const TriKey& k : mScratchFaces)
            if (mStore.mBoundary.at(k).mPatch != mStore.mFeature[a])
                return AdvReject::Class;
    } else if (ka == TetVertexKind::Ridge && mStore.RidgeChain(A, B) != mStore.mFeature[a]) {
        return AdvReject::Class;
    }
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
    if (ka != TetVertexKind::Interior)
        return BoundaryCollapse(A, B);
    return AdvReject::None;
}

void AdvRemesher::CommitCollapse(std::int64_t A, std::int64_t B) {
    const std::size_t a = static_cast<std::size_t>(A);
    if (mStore.Kind(A) != TetVertexKind::Interior)
        CommitBoundaryCollapse(A, B);
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
    const Vec3 mid = detail::vec3_add(pa, detail::vec3_scale(detail::vec3_sub(pb, pa), rT));
    // every child positive and no worse than the floor, with the new vertex at rQ
    auto children = [&](const Vec3& rQ) {
        double old_min = std::numeric_limits<double>::infinity();
        double new_min = std::numeric_limits<double>::infinity();
        for (const std::int64_t t : mScratchShared) {
            const TetCorners& c = mStore.mTets[static_cast<std::size_t>(t)];
            old_min = std::min(old_min, TetQuality(c));
            for (const std::int64_t replaced : {B, A}) {
                Vec3 q[4];
                for (int k = 0; k < 4; ++k)
                    q[k] = c[k] == replaced ? rQ : mStore.Point(c[k]);
                if (detail::tet_volume6(q[0], q[1], q[2], q[3]) <= mEpsVol6 ||
                    !TetRemeshStore::Positive(q[0], q[1], q[2], q[3]))
                    return AdvReject::Inversion;
                new_min = std::min(new_min, Quality(q[0], q[1], q[2], q[3]));
            }
        }
        if (new_min < std::min(old_min, mQualityMin))
            return AdvReject::Quality;
        return AdvReject::None;
    };
    // a boundary edge: its midpoint on the chain, or projected onto its patch,
    // when the children, the normals and the envelope allow; else the midpoint
    if (!mScratchFaces.empty()) {
        Vec3 proj{0.0, 0.0, 0.0};
        bool have = false;
        const std::int32_t chain = mStore.RidgeChain(A, B);
        if (chain >= 0) {
            const double sa = mStore.ChainParam(A, chain, mStore.ChainParam(B, chain, 0.0));
            const double sb = mStore.ChainParam(B, chain, sa);
            proj = mStore.mChains[static_cast<std::size_t>(chain)].Eval(sa + (sb - sa) * rT);
            have = true;
        } else {
            const std::int32_t patch = mStore.mBoundary.at(mScratchFaces.front()).mPatch;
            bool one_patch = true;
            for (const TriKey& k : mScratchFaces)
                one_patch = one_patch && mStore.mBoundary.at(k).mPatch == patch;
            have = one_patch && ProjectToPatch(mid, patch, proj);
        }
        if (have && children(proj) == AdvReject::None) {
            std::vector<std::array<Vec3, 3>> faces;
            std::vector<std::int32_t> patches;
            bool ok = true;
            for (const TriKey& k : mScratchFaces) {
                const auto o = mStore.OrientedFace(k);
                std::array<Vec3, 3> old{};
                for (std::size_t i = 0; i < 3; ++i)
                    old[i] = mStore.Point(o[i]);
                const Vec3 n0 = detail::vec3_cross(detail::vec3_sub(old[1], old[0]),
                                                   detail::vec3_sub(old[2], old[0]));
                for (const std::int64_t replaced : {B, A}) {
                    std::array<Vec3, 3> f = old;
                    for (std::size_t i = 0; i < 3; ++i)
                        if (o[i] == replaced)
                            f[i] = proj;
                    const Vec3 n1 = detail::vec3_cross(detail::vec3_sub(f[1], f[0]),
                                                       detail::vec3_sub(f[2], f[0]));
                    const double l0 = detail::vec3_norm(n0), l1 = detail::vec3_norm(n1);
                    ok = ok && l1 > 0.0 && detail::vec3_dot(n0, n1) >= mCosNormal * l0 * l1;
                    faces.push_back(f);
                    patches.push_back(mStore.mBoundary.at(k).mPatch);
                }
            }
            std::vector<std::int64_t> around{A, B};
            for (const TriKey& k : mScratchFaces)
                around.insert(around.end(), k.begin(), k.end());
            std::sort(around.begin(), around.end());
            around.erase(std::unique(around.begin(), around.end()), around.end());
            if (ok && InEnvelope(faces, patches) &&
                ReverseCovered(FaceGeometry(mScratchFaces), faces, around, mScratchFaces)) {
                rP = proj;
                return AdvReject::None;
            }
        }
    }
    rP = mid;
    return children(mid);
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

// --- the boundary ------------------------------------------------------------------

void AdvRemesher::BuildReference() {
    std::vector<TriKey> keys;
    keys.reserve(mStore.mBoundary.size());
    for (const auto& kv : mStore.mBoundary)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    mCosNormal = std::cos(mrOptions.mMaxNormalDeviation * 3.14159265358979323846 / 180.0);
    if (keys.empty())
        return;
    const std::size_t n = mStore.NumVertices();
    mRefSoup.mPoints.resize(n);
    for (std::size_t v = 0; v < n; ++v)
        mRefSoup.mPoints[v] = mStore.Point(static_cast<std::int64_t>(v));
    for (std::size_t f = 0; f < keys.size(); ++f) {
        const auto o = mStore.OrientedFace(keys[f]);
        mRefTriangles.push_back(o);
        for (const std::int64_t v : o)
            mRefSoup.mCorners.push_back(mStore.Point(v));
        mRefSoup.mVertices.push_back(o);
        mRefSoup.mSourceCell.push_back(static_cast<std::int64_t>(f));
        mRefPatch.push_back(mStore.mBoundary.at(keys[f]).mPatch);
    }
    SurfaceDistanceOptions so;
    so.mWatertightCheck = SdfWatertightCheck::Off;
    mRefQuery = std::make_unique<detail::DistanceQuery>(detail::build_distance_query(mRefSoup, so));

    // reverse samples: every used vertex, edge midpoint and centroid
    std::vector<std::uint8_t> used(n, 0);
    std::vector<std::uint64_t> edges;
    double edge_sum = 0.0;
    for (const auto& t : mRefTriangles) {
        Vec3 c{0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < 3; ++i) {
            used[static_cast<std::size_t>(t[i])] = 1;
            c = detail::vec3_add(c, mStore.Point(t[i]));
            edges.push_back(detail::edge_key64(t[i], t[(i + 1) % 3]));
        }
        mRefSamples.push_back(detail::vec3_scale(c, 1.0 / 3.0));
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    for (const std::uint64_t e : edges) {
        const Vec3 pa = mStore.Point(static_cast<std::int64_t>(e >> 32));
        const Vec3 pb = mStore.Point(static_cast<std::int64_t>(e & 0xFFFFFFFFull));
        mRefSamples.push_back(detail::vec3_scale(detail::vec3_add(pa, pb), 0.5));
        edge_sum += detail::vec3_norm(detail::vec3_sub(pa, pb));
    }
    for (std::size_t v = 0; v < n; ++v)
        if (used[v])
            mRefSamples.push_back(mStore.Point(static_cast<std::int64_t>(v)));
    const double cell =
        edges.empty() ? mDiagonal
                      : std::max(edge_sum / static_cast<double>(edges.size()), 1e-12 * mDiagonal);
    mRefGrid = std::make_unique<detail::SpatialGrid>(cell);
    for (std::size_t i = 0; i < mRefSamples.size(); ++i)
        mRefGrid->Insert(mRefGrid->KeyOf(mRefSamples[i].data()), static_cast<std::int64_t>(i));
}

std::vector<std::array<Vec3, 3>> AdvRemesher::FaceGeometry(const std::vector<TriKey>& rKeys) const {
    std::vector<std::array<Vec3, 3>> out;
    out.reserve(rKeys.size());
    for (const TriKey& k : rKeys)
        out.push_back({mStore.Point(k[0]), mStore.Point(k[1]), mStore.Point(k[2])});
    return out;
}

bool AdvRemesher::ReverseCovered(const std::vector<std::array<Vec3, 3>>& rOld,
                                 const std::vector<std::array<Vec3, 3>>& rNew,
                                 const std::vector<std::int64_t>& rAround,
                                 const std::vector<TriKey>& rReplaced) {
    if (!mRefGrid || rOld.empty())
        return true;
    Vec3 lo = rOld[0][0], hi = rOld[0][0];
    for (const auto& f : rOld)
        for (const Vec3& p : f)
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], p[d]);
                hi[d] = std::max(hi[d], p[d]);
            }
    for (std::size_t d = 0; d < 3; ++d) {
        lo[d] -= mHausd;
        hi[d] += mHausd;
    }
    auto dist2 = [](const Vec3& rP, const std::vector<std::array<Vec3, 3>>& rFaces) {
        double best = std::numeric_limits<double>::infinity();
        for (const auto& f : rFaces)
            best =
                std::min(best, detail::closest_point_on_triangle(rP, f[0], f[1], f[2]).mDistanceSq);
        return best;
    };
    // the unchanged faces around: boundary faces of the vertices around, minus
    // the ones the edit replaces
    std::vector<TriKey> around, faces;
    for (const std::int64_t v : rAround) {
        mStore.BoundaryFacesOf(v, faces);
        for (const TriKey& k : faces)
            if (!std::binary_search(rReplaced.begin(), rReplaced.end(), k))
                around.push_back(k);
    }
    std::sort(around.begin(), around.end());
    around.erase(std::unique(around.begin(), around.end()), around.end());
    std::vector<std::array<Vec3, 3>> keep = rNew;
    const auto unchanged = FaceGeometry(around);
    keep.insert(keep.end(), unchanged.begin(), unchanged.end());
    const double h2 = mHausd * mHausd, margin2 = 0.81 * h2;
    bool ok = true;
    std::vector<std::int64_t> seen;
    mRefGrid->ForEachInBox(mRefGrid->KeyOf(lo.data()), mRefGrid->KeyOf(hi.data()),
                           [&](const std::vector<std::int64_t>& rIds) {
                               for (const std::int64_t i : rIds) {
                                   if (!ok)
                                       return;
                                   const Vec3& p = mRefSamples[static_cast<std::size_t>(i)];
                                   bool inside = true;
                                   for (std::size_t d = 0; d < 3; ++d)
                                       inside = inside && p[d] >= lo[d] && p[d] <= hi[d];
                                   if (!inside)
                                       continue;
                                   const double before = dist2(p, rOld);
                                   if (before > h2)
                                       continue;
                                   if (dist2(p, keep) > std::max(margin2, before))
                                       ok = false;
                               }
                           });
    return ok;
}

bool AdvRemesher::ProjectToPatch(const Vec3& rP, std::int32_t Patch, Vec3& rOut) const {
    if (!mRefQuery)
        return false;
    const auto hits = detail::query_surface_projections_filtered(
        *mRefQuery, {rP},
        [&](std::size_t, std::int64_t T) {
            return mRefPatch[static_cast<std::size_t>(T)] == Patch;
        },
        std::numeric_limits<double>::infinity());
    if (!hits[0].mFound)
        return false;
    rOut = hits[0].mPoint;
    return true;
}

bool AdvRemesher::InEnvelope(const std::vector<std::array<Vec3, 3>>& rTris,
                             const std::vector<std::int32_t>& rPatches) const {
    if (!mRefQuery || rTris.empty())
        return true;
    // a barycentric grid per face, its spacing at most hausd (2 to 6
    // subdivisions), every sample within 0.9 hausd of the face's patch: the
    // margin covers the face between the samples
    std::vector<Vec3> pts;
    std::vector<std::int32_t> patch_of;
    for (std::size_t f = 0; f < rTris.size(); ++f) {
        const auto& t = rTris[f];
        double lmax = 0.0;
        for (std::size_t i = 0; i < 3; ++i)
            lmax = std::max(lmax, detail::vec3_norm(detail::vec3_sub(t[(i + 1) % 3], t[i])));
        const int n = std::clamp(static_cast<int>(std::ceil(lmax / mHausd)), 2, 6);
        for (int i = 0; i <= n; ++i)
            for (int j = 0; i + j <= n; ++j) {
                const double u = static_cast<double>(i) / n, v = static_cast<double>(j) / n;
                const double w = 1.0 - u - v;
                pts.push_back(detail::vec3_add(
                    detail::vec3_add(detail::vec3_scale(t[0], w), detail::vec3_scale(t[1], u)),
                    detail::vec3_scale(t[2], v)));
                patch_of.push_back(rPatches[f]);
            }
    }
    const auto hits = detail::query_surface_projections_filtered(
        *mRefQuery, pts,
        [&](std::size_t Q, std::int64_t T) {
            return mRefPatch[static_cast<std::size_t>(T)] == patch_of[Q];
        },
        0.9 * mHausd);
    for (const auto& h : hits)
        if (!h.mFound)
            return false;
    return true;
}

bool AdvRemesher::MovedFaces(std::int64_t V, std::int64_t Onto, const Vec3& rP,
                             std::vector<std::array<Vec3, 3>>& rNew,
                             std::vector<std::int32_t>& rPatches) {
    rNew.clear();
    rPatches.clear();
    mStore.BoundaryFacesOf(V, mScratchTriA);
    for (const TriKey& k : mScratchTriA) {
        if (Onto >= 0 && (k[0] == Onto || k[1] == Onto || k[2] == Onto))
            continue;  // vanishes with the collapse
        const auto o = mStore.OrientedFace(k);
        std::array<Vec3, 3> before{}, after{};
        for (std::size_t i = 0; i < 3; ++i) {
            before[i] = mStore.Point(o[i]);
            after[i] = o[i] == V ? rP : before[i];
        }
        const Vec3 n0 = detail::vec3_cross(detail::vec3_sub(before[1], before[0]),
                                           detail::vec3_sub(before[2], before[0]));
        const Vec3 n1 = detail::vec3_cross(detail::vec3_sub(after[1], after[0]),
                                           detail::vec3_sub(after[2], after[0]));
        const double l0 = detail::vec3_norm(n0), l1 = detail::vec3_norm(n1);
        if (!(l1 > 0.0))
            return false;
        if (l0 > 0.0 && detail::vec3_dot(n0, n1) < mCosNormal * l0 * l1)
            return false;
        rNew.push_back(after);
        rPatches.push_back(mStore.mBoundary.at(k).mPatch);
    }
    return true;
}

// The boundary half of a collapse of a boundary vertex a onto b: the boundary
// link condition (the vertices sharing a boundary face with both a and b are
// exactly the far corners of the boundary faces through ab, and no two
// boundary faces become one), the normals, and the envelope both ways.
AdvReject AdvRemesher::BoundaryCollapse(std::int64_t A, std::int64_t B) {
    mStore.BoundaryFacesOf(A, mScratchTriA);
    mStore.BoundaryFacesOf(B, mScratchTriB);
    auto ring = [](const std::vector<TriKey>& rFaces, std::int64_t V) {
        std::vector<std::int64_t> r;
        for (const TriKey& k : rFaces)
            for (const std::int64_t x : k)
                if (x != V)
                    r.push_back(x);
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        return r;
    };
    const std::vector<std::int64_t> ra = ring(mScratchTriA, A), rb = ring(mScratchTriB, B);
    std::vector<std::int64_t> common, opposite;
    std::set_intersection(ra.begin(), ra.end(), rb.begin(), rb.end(), std::back_inserter(common));
    for (const TriKey& k : mScratchTriA)
        if (k[0] == B || k[1] == B || k[2] == B)
            for (const std::int64_t x : k)
                if (x != A && x != B)
                    opposite.push_back(x);
    std::sort(opposite.begin(), opposite.end());
    opposite.erase(std::unique(opposite.begin(), opposite.end()), opposite.end());
    common.erase(std::remove(common.begin(), common.end(), A), common.end());
    common.erase(std::remove(common.begin(), common.end(), B), common.end());
    if (common != opposite)
        return AdvReject::Link;
    for (const TriKey& k : mScratchTriA) {
        if (k[0] == B || k[1] == B || k[2] == B)
            continue;
        std::int64_t o[2];
        int i = 0;
        for (const std::int64_t x : k)
            if (x != A)
                o[i++] = x;
        if (mStore.mBoundary.count(detail::tri_key(B, o[0], o[1])))
            return AdvReject::Link;
    }
    std::vector<std::array<Vec3, 3>> faces;
    std::vector<std::int32_t> patches;
    if (!MovedFaces(A, B, mStore.Point(B), faces, patches))
        return AdvReject::Normal;
    if (!InEnvelope(faces, patches))
        return AdvReject::Hausdorff;
    // the other way: the reference the old faces covered stays covered
    mStore.BoundaryFacesOf(A, mScratchTriA);
    if (!ReverseCovered(FaceGeometry(mScratchTriA), faces, ra, mScratchTriA))
        return AdvReject::Hausdorff;
    return AdvReject::None;
}

void AdvRemesher::CommitBoundaryCollapse(std::int64_t A, std::int64_t B) {
    std::vector<TriKey> faces;
    mStore.BoundaryFacesOf(A, faces);
    for (const TriKey& k : faces) {
        const detail::TetBoundaryFace face = mStore.mBoundary.at(k);
        mStore.mBoundary.erase(k);
        if (k[0] == B || k[1] == B || k[2] == B)
            continue;
        std::int64_t o[2];
        int i = 0;
        for (const std::int64_t x : k)
            if (x != A)
                o[i++] = x;
        mStore.mBoundary.emplace(detail::tri_key(B, o[0], o[1]), face);
    }
    std::vector<std::int64_t> ring;
    mStore.VertexRing(A, ring);
    for (const std::int64_t x : ring) {
        const std::int32_t chain = mStore.RidgeChain(A, x);
        if (chain < 0)
            continue;
        mStore.mRidgeEdges.erase(detail::edge_key64(A, x));
        if (x != B)
            mStore.mRidgeEdges[detail::edge_key64(B, x)] = chain;
    }
}

double AdvRemesher::StarMinQuality(std::int64_t V) const {
    double m = std::numeric_limits<double>::infinity();
    for (const std::int64_t t : mStore.mVertexTets[static_cast<std::size_t>(V)])
        m = std::min(m, TetQuality(mStore.mTets[static_cast<std::size_t>(t)]));
    return m;
}

bool AdvRemesher::StarValid(std::int64_t V, const Vec3& rP, double Floor) const {
    for (const std::int64_t t : mStore.mVertexTets[static_cast<std::size_t>(V)]) {
        const TetCorners& c = mStore.mTets[static_cast<std::size_t>(t)];
        Vec3 q[4];
        for (int k = 0; k < 4; ++k)
            q[k] = c[k] == V ? rP : mStore.Point(c[k]);
        if (detail::tet_volume6(q[0], q[1], q[2], q[3]) <= mEpsVol6 ||
            !TetRemeshStore::Positive(q[0], q[1], q[2], q[3]))
            return false;
        if (Quality(q[0], q[1], q[2], q[3]) < Floor)
            return false;
    }
    return true;
}

// Surface vertices slide half way to the mean of their boundary neighbours
// along the tangent plane and back onto their patch of the reference skin;
// ridge vertices slide half way to the middle of their two chain neighbours,
// along the chain. A move is kept when the star's least quality does not
// drop, every normal turns by at most the bound and the faces stay within
// hausd of their patch.
std::int64_t AdvRemesher::BoundarySmoothPass() {
    std::int64_t moved = 0;
    std::vector<TriKey> faces;
    std::vector<std::array<Vec3, 3>> new_faces;
    std::vector<std::int32_t> patches;
    std::vector<std::int64_t> ring;
    for (std::size_t v = 0; v < mStore.NumVertices(); ++v) {
        const std::int64_t V = static_cast<std::int64_t>(v);
        const TetVertexKind kind = mStore.Kind(V);
        if (!mStore.mVertexAlive[v] || mStore.mFrozen[v] ||
            (kind != TetVertexKind::Surface && kind != TetVertexKind::Ridge))
            continue;
        const Vec3 p0 = mStore.Point(V);
        Vec3 target{0.0, 0.0, 0.0};
        double param = mStore.mChainParam[v];
        if (kind == TetVertexKind::Surface) {
            mStore.BoundaryFacesOf(V, faces);
            Vec3 normal{0.0, 0.0, 0.0}, mean{0.0, 0.0, 0.0};
            std::vector<std::int64_t> nb;
            for (const TriKey& k : faces) {
                const auto o = mStore.OrientedFace(k);
                normal = detail::vec3_add(
                    normal,
                    detail::vec3_cross(detail::vec3_sub(mStore.Point(o[1]), mStore.Point(o[0])),
                                       detail::vec3_sub(mStore.Point(o[2]), mStore.Point(o[0]))));
                for (const std::int64_t x : k)
                    if (x != V)
                        nb.push_back(x);
            }
            std::sort(nb.begin(), nb.end());
            nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
            if (nb.empty())
                continue;
            for (const std::int64_t x : nb)
                mean = detail::vec3_add(mean, mStore.Point(x));
            mean = detail::vec3_scale(mean, 1.0 / static_cast<double>(nb.size()));
            const Vec3 n = detail::vec3_normalize(normal);
            Vec3 d = detail::vec3_scale(detail::vec3_sub(mean, p0), 0.5);
            d = detail::vec3_sub(d, detail::vec3_scale(n, detail::vec3_dot(d, n)));
            if (!ProjectToPatch(detail::vec3_add(p0, d), mStore.mFeature[v], target))
                continue;
        } else {
            const std::int32_t chain = mStore.mFeature[v];
            mStore.VertexRing(V, ring);
            std::vector<double> s;
            for (const std::int64_t x : ring)
                if (mStore.RidgeChain(V, x) == chain)
                    s.push_back(mStore.ChainParam(x, chain, param));
            if (s.size() != 2)
                continue;
            param = param + 0.5 * (0.5 * (s[0] + s[1]) - param);
            target = mStore.mChains[static_cast<std::size_t>(chain)].Eval(param);
        }
        if (detail::vec3_norm_sq(detail::vec3_sub(target, p0)) <= 0.0)
            continue;
        if (!StarValid(V, target, StarMinQuality(V)))
            continue;
        if (!MovedFaces(V, -1, target, new_faces, patches) || !InEnvelope(new_faces, patches))
            continue;
        mStore.BoundaryFacesOf(V, faces);
        std::vector<std::int64_t> around;
        for (const TriKey& k : faces)
            around.insert(around.end(), k.begin(), k.end());
        std::sort(around.begin(), around.end());
        around.erase(std::unique(around.begin(), around.end()), around.end());
        if (!ReverseCovered(FaceGeometry(faces), new_faces, around, faces))
            continue;
        mStore.SetPoint(V, target);
        if (kind == TetVertexKind::Ridge)
            mStore.mChainParam[v] = param;
        ++mVersion[v];
        ++moved;
    }
    return moved;
}

double AdvRemesher::BoundaryDeviation() const {
    if (mRefTriangles.empty())
        return 0.0;
    auto skin = [](const std::vector<double>& rXyz, std::size_t N,
                   const std::vector<std::array<std::int64_t, 3>>& rTris) {
        Mesh m;
        NDArray pts = NDArray::Uninit(DType::Float64, {N, 3});
        std::copy(rXyz.begin(), rXyz.begin() + static_cast<std::ptrdiff_t>(N * 3),
                  pts.As<double>());
        m.AssignPoints(std::move(pts));
        NDArray conn = NDArray::Uninit(DType::Int64, {rTris.size(), 3});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rTris.size(); ++i)
            for (std::size_t k = 0; k < 3; ++k)
                c[i * 3 + k] = rTris[i][k];
        m.AddCellBlock(cell_type_name(CellType::Triangle), std::move(conn));
        return m;
    };
    std::vector<double> ref(mRefSoup.mPoints.size() * 3);
    for (std::size_t v = 0; v < mRefSoup.mPoints.size(); ++v)
        for (std::size_t d = 0; d < 3; ++d)
            ref[v * 3 + d] = mRefSoup.mPoints[v][d];
    std::vector<TriKey> keys;
    for (const auto& kv : mStore.mBoundary)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    std::vector<std::array<std::int64_t, 3>> now;
    for (const TriKey& k : keys)
        now.push_back(mStore.OrientedFace(k));
    HausdorffOptions ho;
    ho.mFaceSamples = 4;
    const HausdorffResult h =
        hausdorff_distance(skin(mStore.mXyz, mStore.NumVertices(), now),
                           skin(ref, mRefSoup.mPoints.size(), mRefTriangles), ho);
    return h.mDistance;
}

// --- smoothing -------------------------------------------------------------------

std::vector<std::uint8_t> AdvRemesher::InteriorPins(std::int32_t Class) const {
    std::vector<std::uint8_t> pinned(mStore.NumVertices(), 1);
    for (std::size_t v = 0; v < mStore.NumVertices(); ++v) {
        if (!mStore.mVertexAlive[v] || mStore.mFrozen[v] ||
            mStore.Kind(static_cast<std::int64_t>(v)) != TetVertexKind::Interior ||
            mStore.mVertexTets[v].empty())
            continue;
        if (mStore.mClass[static_cast<std::size_t>(mStore.mVertexTets[v].front())] == Class)
            pinned[v] = 0;
    }
    return pinned;
}

// The interior vertices of each class: optimize_volume's guarded ODT
// relocation (each moved vertex's star keeps its least dihedral angle) and,
// for the tetrahedra under 40 degrees, its worst-cell smoothing.
std::int64_t AdvRemesher::SmoothInteriorPass(bool Odt) {
    std::vector<std::vector<detail::TetFlipTet>> by_class;
    std::vector<detail::TetFlipTet> fixed;
    std::vector<std::int32_t> fixed_class;
    TetsByClass(by_class, fixed, fixed_class);
    SmoothOptions so;
    so.mMethod = SmoothMethod::Odt;
    so.mIterations = 1;
    so.mGuardInversion = true;
    const std::vector<double> before = mStore.mXyz;
    for (std::size_t c = 0; c < by_class.size(); ++c) {
        if (by_class[c].empty())
            continue;
        const std::vector<std::uint8_t> pinned = InteriorPins(static_cast<std::int32_t>(c));
        if (std::find(pinned.begin(), pinned.end(), 0) == pinned.end())
            continue;
        if (Odt)
            detail::tet_flip_guarded_odt(mStore.mXyz, by_class[c], pinned, so, mEpsVol6);
        detail::tet_flip_smooth_worst(mStore.mXyz, by_class[c], pinned, mEpsVol6,
                                      kAdvWorstThresholdDeg);
    }
    // The shared kernels guard dihedral angles, which do not see a
    // tetrahedron's sign: a tetrahedron that was flat can come out inverted.
    // Put back, in ascending order and until none is left, the moved vertices
    // of every tetrahedron the exact orient3d finds inverted, or flat where it
    // was not.
    auto moved_vertex = [&](std::int64_t V) {
        const std::size_t i = static_cast<std::size_t>(V) * 3;
        return mStore.mXyz[i] != before[i] || mStore.mXyz[i + 1] != before[i + 1] ||
               mStore.mXyz[i + 2] != before[i + 2];
    };
    for (bool again = true; again;) {
        again = false;
        for (std::size_t t = 0; t < mStore.NumTets(); ++t) {
            if (!mStore.mTetAlive[t])
                continue;
            const TetCorners& c = mStore.mTets[t];
            bool any_moved = false;
            for (const std::int64_t v : c)
                any_moved = any_moved || moved_vertex(v);
            if (!any_moved)
                continue;
            const int sign = TetRemeshStore::OrientSign(mStore.Point(c[0]), mStore.Point(c[1]),
                                                        mStore.Point(c[2]), mStore.Point(c[3]));
            if (sign > 0)
                continue;
            if (sign == 0) {
                TetCorners k = c;
                std::sort(k.begin(), k.end());
                if (std::binary_search(mFlatKeys.begin(), mFlatKeys.end(), k)) {
                    // flat before too? only then may it stay flat
                    const auto pb = [&](std::int64_t V) {
                        const std::size_t i = static_cast<std::size_t>(V) * 3;
                        return Vec3{before[i], before[i + 1], before[i + 2]};
                    };
                    if (TetRemeshStore::OrientSign(pb(c[0]), pb(c[1]), pb(c[2]), pb(c[3])) == 0)
                        continue;
                }
            }
            for (const std::int64_t v : c)
                if (moved_vertex(v))
                    for (std::size_t d = 0; d < 3; ++d)
                        mStore.mXyz[static_cast<std::size_t>(v) * 3 + d] =
                            before[static_cast<std::size_t>(v) * 3 + d];
            again = true;
        }
    }
    std::int64_t moved = 0;
    for (std::size_t v = 0; v < mStore.NumVertices(); ++v)
        moved += moved_vertex(static_cast<std::int64_t>(v));
    if (moved > 0)
        for (std::uint32_t& v : mVersion)
            ++v;
    return moved;
}

// The AMIPS energy of the star of V with V at rP, its least quality, and
// whether every tetrahedron of it stays strictly positive (exact).
double AdvRemesher::StarAmips(std::int64_t V, const Vec3& rP, double& rMinQuality,
                              bool& rValid) const {
    double e = 0.0;
    rMinQuality = std::numeric_limits<double>::infinity();
    rValid = true;
    for (const std::int64_t t : mStore.mVertexTets[static_cast<std::size_t>(V)]) {
        const TetCorners& c = mStore.mTets[static_cast<std::size_t>(t)];
        Vec3 q[4];
        for (int k = 0; k < 4; ++k)
            q[k] = c[k] == V ? rP : mStore.Point(c[k]);
        const double qual = Quality(q[0], q[1], q[2], q[3]);
        if (!(qual > 0.0)) {
            rValid = false;
            return std::numeric_limits<double>::infinity();
        }
        rMinQuality = std::min(rMinQuality, qual);
        e += detail::tet_amips(q[0], q[1], q[2], q[3]);
    }
    return e;
}

// Interior vertices, one at a time in ascending order: a step down the
// gradient of the star's conformal AMIPS energy (central differences), halved
// until the energy drops with every tetrahedron positive and the star's least
// quality no lower (Hu et al. 2020 smooth the same energy).
std::int64_t AdvRemesher::AmipsPass() {
    std::int64_t moved = 0;
    std::vector<std::int64_t> ring;
    for (std::size_t v = 0; v < mStore.NumVertices(); ++v) {
        const std::int64_t V = static_cast<std::int64_t>(v);
        if (!mStore.mVertexAlive[v] || mStore.mFrozen[v] ||
            mStore.Kind(V) != TetVertexKind::Interior)
            continue;
        const Vec3 p0 = mStore.Point(V);
        double q0 = 0.0;
        bool valid = false;
        const double e0 = StarAmips(V, p0, q0, valid);
        if (!valid)
            continue;
        mStore.VertexRing(V, ring);
        double edge = 0.0;
        for (const std::int64_t x : ring)
            edge += detail::vec3_norm(detail::vec3_sub(mStore.Point(x), p0));
        edge /= static_cast<double>(std::max<std::size_t>(ring.size(), 1));
        const double h = 1e-6 * edge;
        Vec3 g{0.0, 0.0, 0.0};
        bool ok = true;
        for (int d = 0; d < 3 && ok; ++d) {
            Vec3 pp = p0, pm = p0;
            pp[static_cast<std::size_t>(d)] += h;
            pm[static_cast<std::size_t>(d)] -= h;
            double qa = 0.0, qb = 0.0;
            bool va = false, vb = false;
            const double ea = StarAmips(V, pp, qa, va), eb = StarAmips(V, pm, qb, vb);
            ok = va && vb;
            g[static_cast<std::size_t>(d)] = (ea - eb) / (2.0 * h);
        }
        const double gn = detail::vec3_norm(g);
        if (!ok || !(gn > 0.0))
            continue;
        double step = 0.1 * edge / gn;
        for (int k = 0; k < 12; ++k, step *= 0.5) {
            const Vec3 p = detail::vec3_sub(p0, detail::vec3_scale(g, step));
            double q = 0.0;
            bool vp = false;
            const double e = StarAmips(V, p, q, vp);
            if (vp && e < e0 && q >= q0) {
                mStore.SetPoint(V, p);
                ++mVersion[v];
                ++moved;
                break;
            }
        }
    }
    return moved;
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
    BuildReference();
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
        if (mrOptions.mEnableSmooth) {
            t0 = AdvClock::now();
            const bool amips = mrOptions.mQualityMetric == AdaptVolumeQuality::Amips;
            std::int64_t n = amips ? AmipsPass() : 0;
            Verify("amips smoothing");
            n += SmoothInteriorPass(!amips);
            Verify("interior smoothing");
            n += BoundarySmoothPass();
            mrResult.mNumSmoothed += n;  // moves do not count towards the 1% stop
            Time("smooth", t0);
            Verify("smooth");
        }
        ++mrResult.mIterations;
        if (static_cast<double>(changed) < 0.01 * static_cast<double>(mStore.NumAliveTets()))
            break;
        if (mrOptions.mTargetMinDihedral > 0.0 &&
            Measure().mMinDihedral >= mrOptions.mTargetMinDihedral)
            break;
    }

    // polish: the worst-cell moves and smoothing until they stop (at most 3 rounds)
    if (mrOptions.mEnableSwap || mrOptions.mEnableSmooth) {
        t0 = AdvClock::now();
        for (int round = 0; round < 3; ++round) {
            std::int64_t n = 0;
            if (mrOptions.mEnableSwap) {
                const std::int64_t k = SwapPass();
                mrResult.mNumSwapped += k;
                n += k;
            }
            if (mrOptions.mEnableSmooth) {
                const std::int64_t k = SmoothInteriorPass(false);
                mrResult.mNumSmoothed += k;
                n += k;
            }
            if (n == 0)
                break;
        }
        Time("polish", t0);
        Verify("polish");
    }

    t0 = AdvClock::now();
    const AdvMeasures after = Measure();
    mrResult.mMinDihedralAfter = after.mMinDihedral;
    mrResult.mMaxDihedralAfter = after.mMaxDihedral;
    mrResult.mMeanQualityAfter = after.mMeanQuality;
    mrResult.mVolumeAfter = after.mVolume;
    mrResult.mMaxBoundaryDeviation = BoundaryDeviation();
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
