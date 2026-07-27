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
#pragma once

/**
 * @file kratos_mesh.hpp
 * @brief The KRATOS mesh backend: `meshioplusplus::KratosMesh`, the uniform
 * mesh API implemented over a Kratos-style `ModelPart`.
 *
 * Selected by `MESHIOPLUSPLUS_MESH_BACKEND=KRATOS` (see `mesh.hpp`). Format
 * readers ingest into a canonical staging structure (a `NativeMesh` — same
 * canonical Float64/Int64 storage), and the `ModelPart` is **materialized
 * lazily** on the first `GetModelPart()` call:
 *
 *  - Nodes get Ids `index + 1` (z = 0-padded for 2-D points).
 *  - Cell blocks whose topological dimension equals the mesh's maximum
 *    become **Elements**; lower-dimension blocks become **Conditions** (the
 *    Kratos convention, matching `src/python/meshioplusplus/mdpa/_mdpa.py`), each
 *    kind Id-numbered 1..N in block order, with default Kratos names from
 *    `kratos_names.hpp`.
 *  - `point_data` becomes nodal data; `cell_data` becomes elemental /
 *    conditional data (concatenated per kind, entity order); `field_data`
 *    stays on the staging mesh (Kratos has no equivalent).
 *  - **Integer tag arrays** under well-known names (`gmsh:physical`,
 *    `su2:tag`, `medit:ref`, `cell_tags`, ...) automatically become named
 *    SubModelParts (`gmsh_physical_1`, ...) containing the tagged entities
 *    and their nodes — disable with `SetBuildSubModelPartsFromTags(false)`.
 *    The tag arrays remain as elemental/conditional data either way, so
 *    writer round-trips are unaffected.
 *  - **Ragged blocks** (polygon/polyhedron) have no ModelPart geometry;
 *    they stay in staging pass-through (round-trips keep working) and do
 *    not become entities.
 *
 * Writer accessors serve from the staging structure, so a read -> write
 * round-trip never pays for (or builds) the ModelPart at all, and output
 * bytes match the NATIVE backend exactly. After mutating the ModelPart
 * directly, call `InvalidateBlocks()`: the staging is then rebuilt from the
 * ModelPart on the next accessor use (consecutive same-type Elements are
 * grouped into blocks, then Conditions). Since v9.0.0 the rebuild is
 * lossless for everything the ModelPart can express — ragged pass-through
 * blocks are carried through at their original positions (they never became
 * entities, so a ModelPart edit cannot have touched them) and SubModelParts
 * are read back as named `Cell`/`Point` regions. Two things still cannot
 * survive, because a SubModelPart has nowhere to put them: a region's
 * `mDim`/`mTag` (recovered from a staged region of the same name when there is
 * one, else `-1`) and SubModelPart *nesting*, regions being flat.
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/backends/model_part.hpp"
#include "meshioplusplus/backends/native_mesh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

/**
 * @brief The KRATOS mesh backend (aliased to `meshioplusplus::Mesh` when
 * `MESHIOPLUSPLUS_MESH_BACKEND_KRATOS` is defined). See the file-level
 * comment for the staging/materialization design.
 */
class KratosMesh {
public:
    /**
     * @brief The cell_data key holding each cell's Kratos properties id.
     *
     * The same key the MDPA reader/writer uses, so a `.mdpa` round trip through
     * this backend keeps its material assignment instead of collapsing to 0.
     */
    static constexpr const char* kKratosPropertyTagKey = "gmsh:physical";

    // --- KRATOS-only extras: per-block Kratos entity names -----------------
    //
    // Not part of the uniform mesh API (no other backend has a ModelPart to
    // spell names for) -- the same "fast-consumer surface" shape as NATIVE's
    // PointsData()/ConnSpan(). A block with no name set derives one from its
    // cell type, which is what every caller got before names existed.

    /** @brief Set block @p Block's Kratos entity name (e.g. from `MdpaInfo`). */
    void SetBlockEntityName(std::size_t Block, const std::string& rName) {
        if (Block >= mStage.NumCellBlocks())
            throw std::out_of_range("meshio++ KratosMesh: block index out of range");
        if (mBlockEntityNames.size() < mStage.NumCellBlocks())
            mBlockEntityNames.resize(mStage.NumCellBlocks());
        mBlockEntityNames[Block] = rName;
        InvalidateMaterialization();
    }
    /** @brief Block @p Block's Kratos entity name, or `""` when unset. */
    const std::string& BlockEntityName(std::size_t Block) const {
        static const std::string empty;
        return Block < mBlockEntityNames.size() ? mBlockEntityNames[Block] : empty;
    }

    /** @brief Cell-data names treated as entity tags for automatic SubModelParts. */
    static const std::vector<std::string>& KnownTagKeys() {
        static const std::vector<std::string> keys = {
            "cell_tags", "gmsh:physical", "su2:tag", "medit:ref",  "avsucd:material", "freefem:ref",
            "mfm:ref",   "netgen:index",  "pf3:ref", "tetgen:ref", "ugrid:ref",       "unv:pid",
        };
        return keys;
    }

    // --- uniform API: reader-side ingestion (forwarded to staging) ---------

    void AssignPoints(NDArray points) {
        ResetModelPartOnly();
        mStage.AssignPoints(std::move(points));
    }
    void AddCellBlock(std::string type, NDArray conn) {
        ResetModelPartOnly();
        mStage.AddCellBlock(std::move(type), std::move(conn));
    }
    void AddPolygonBlock(std::string type, std::vector<std::vector<std::int64_t>> rows) {
        ResetModelPartOnly();
        mStage.AddPolygonBlock(std::move(type), std::move(rows));
    }
    void AddPolyhedronBlock(std::string type,
                            std::vector<std::vector<std::vector<std::int64_t>>> cells) {
        ResetModelPartOnly();
        mStage.AddPolyhedronBlock(std::move(type), std::move(cells));
    }
    void AddPointData(std::string name, NDArray data) {
        ResetModelPartOnly();
        mStage.AddPointData(std::move(name), std::move(data));
    }
    void AddCellData(std::string name, std::vector<NDArray> blocks) {
        ResetModelPartOnly();
        mStage.AddCellData(std::move(name), std::move(blocks));
    }
    void AppendCellData(const std::string& rName, NDArray block) {
        ResetModelPartOnly();
        mStage.AppendCellData(rName, std::move(block));
    }
    void AddFieldData(std::string name, NDArray data) {
        ResetModelPartOnly();
        mStage.AddFieldData(std::move(name), std::move(data));
    }
    void AddRegion(meshioplusplus::Region region) {
        ResetModelPartOnly();
        mStage.AddRegion(std::move(region));
    }
    void AddPropertySet(PropertySet propertySet) {
        ResetModelPartOnly();
        mStage.AddPropertySet(std::move(propertySet));
    }

    // --- uniform API: writer-side accessors (staging, stale-synced) --------

    using CellView = NativeMesh::CellView;

    std::size_t NumPoints() const { return Stage().NumPoints(); }
    std::size_t PointDim() const { return Stage().PointDim(); }
    const NDArray& Points() const { return Stage().Points(); }
    std::size_t NumCellBlocks() const { return Stage().NumCellBlocks(); }
    CellView Cells(std::size_t i) const { return Stage().Cells(i); }
    detail::CellBlockRange<KratosMesh> CellRange() const {
        return detail::CellBlockRange<KratosMesh>(*this);
    }

    std::vector<std::string> PointDataNames() const { return Stage().PointDataNames(); }
    std::size_t NumPointData() const { return Stage().NumPointData(); }
    bool HasPointData(const std::string& rName) const { return Stage().HasPointData(rName); }
    const NDArray& PointData(const std::string& rName) const { return Stage().PointData(rName); }

    std::vector<std::string> CellDataNames() const { return Stage().CellDataNames(); }
    std::size_t NumCellData() const { return Stage().NumCellData(); }
    bool HasCellData(const std::string& rName) const { return Stage().HasCellData(rName); }
    const NDArray& CellData(const std::string& rName, std::size_t block) const {
        return Stage().CellData(rName, block);
    }
    std::size_t CellDataNumBlocks(const std::string& rName) const {
        return Stage().CellDataNumBlocks(rName);
    }

    std::vector<std::string> FieldDataNames() const { return Stage().FieldDataNames(); }
    std::size_t NumFieldData() const { return Stage().NumFieldData(); }
    bool HasFieldData(const std::string& rName) const { return Stage().HasFieldData(rName); }
    const NDArray& FieldData(const std::string& rName) const { return Stage().FieldData(rName); }

    // --- named regions (region.hpp) ----------------------------------------
    // NOTE: `Region(i)` hides the type name `meshioplusplus::Region` for the
    // rest of this class body — every declaration below must qualify it.

    /// Sentinel returned by `FindRegion` when no region matches.
    static constexpr std::size_t npos = detail::RegionList::npos;

    std::size_t NumPropertySets() const { return Stage().NumPropertySets(); }
    const PropertySet& GetPropertySet(std::size_t Index) const {
        return Stage().GetPropertySet(Index);
    }
    bool HasPropertySet(std::int64_t Id) const { return Stage().HasPropertySet(Id); }
    std::size_t FindPropertySet(std::int64_t Id) const { return Stage().FindPropertySet(Id); }

    std::size_t NumRegions() const { return Stage().NumRegions(); }
    const meshioplusplus::Region& Region(std::size_t i) const { return Stage().Region(i); }
    std::vector<std::string> RegionNames() const { return Stage().RegionNames(); }
    bool HasRegion(const std::string& rName) const { return Stage().HasRegion(rName); }
    bool HasRegion(const std::string& rName, RegionKind kind) const {
        return Stage().HasRegion(rName, kind);
    }
    std::size_t FindRegion(const std::string& rName, RegionKind kind) const {
        return Stage().FindRegion(rName, kind);
    }

    // --- KRATOS-specific surface -------------------------------------------

    /**
     * @brief The ModelPart view of this mesh, materialized on first call.
     * @return The root `ModelPart` (named "Main").
     */
    ModelPart& GetModelPart() { return EnsureMaterialized(); }

    /**
     * @brief `const` overload, for wrappers that take a `const Mesh&`.
     *
     * Materialization is lazy, so this DOES do work on first call (and is
     * therefore not thread-safe against a concurrent first call) -- the
     * ModelPart is a cache, exactly like the staging `mStage` that every other
     * const accessor already builds on demand, which is why the members are
     * `mutable`. Logical constness holds: the mesh's value is unchanged.
     */
    const ModelPart& GetModelPart() const { return EnsureMaterialized(); }

    /** @brief Whether `GetModelPart()` has materialized the ModelPart yet. */
    bool IsMaterialized() const { return mMaterialized; }

    /**
     * @brief Declare that the ModelPart was mutated directly: the block/point
     * staging is rebuilt from the ModelPart on the next accessor use.
     *
     * The rebuild groups consecutive same-type Elements, then Conditions, so a
     * mesh whose entity structure is unchanged comes back with its original
     * block layout. Ragged pass-through blocks (which never became entities and
     * so cannot have been edited) are carried through at their original
     * positions, and SubModelParts are read back as named `Cell`/`Point`
     * regions -- neither is lost, as they were before v9.0.0.
     *
     * What still cannot survive: a SubModelPart's *nesting* (regions are flat),
     * and a region's `mDim`/`mTag`, which a SubModelPart has nowhere to store
     * -- a region reconstructed from one comes back with `-1` for both unless a
     * staged region of the same name and kind supplied them.
     */
    void InvalidateBlocks() {
        if (mMaterialized)
            mStale = true;
    }

    /**
     * @brief Enable/disable automatic tag -> SubModelPart creation at
     * materialization (default: enabled). Call before `GetModelPart()`.
     */
    void SetBuildSubModelPartsFromTags(bool enable) {
        mTagsToSubModelParts = enable;
        InvalidateMaterialization();  // re-materialize with the new setting
    }
    bool BuildSubModelPartsFromTags() const { return mTagsToSubModelParts; }

    /**
     * @brief Restrict the automatic tag pass to these `cell_data` keys.
     *
     * `SetBuildSubModelPartsFromTags` is all-or-nothing, which is too coarse
     * when only one key is unwanted. The motivating case is `.mdpa`: a Kratos
     * properties id is read as a `gmsh:physical` cell tag, so the tag pass
     * synthesizes a `gmsh_physical_<id>` SubModelPart beside the file's real
     * ones -- material assignment surfacing as a group. For a genuine gmsh file
     * that same inference is wanted, so the key cannot simply be dropped from
     * #KnownTagKeys; the consumer has to say which it means.
     *
     * @param keys The keys to consider. **Empty restores the default**, i.e.
     *        every key in #KnownTagKeys.
     */
    void SetTagSubModelPartKeys(std::vector<std::string> keys) {
        mTagKeys = std::move(keys);
        InvalidateMaterialization();
    }

    /**
     * @brief Exclude one key from the automatic tag pass. Repeatable, additive.
     *
     * Expresses "everything except this" without re-listing #KnownTagKeys,
     * which an inclusion list alone would force -- and which would silently
     * freeze the caller against any key added to that table later.
     */
    void ExcludeTagSubModelPartKey(const std::string& rKey) {
        if (std::find(mExcludedTagKeys.begin(), mExcludedTagKeys.end(), rKey) ==
            mExcludedTagKeys.end())
            mExcludedTagKeys.push_back(rKey);
        InvalidateMaterialization();
    }

    /**
     * @brief The keys the tag pass will actually consider.
     * @return The effective set, in #KnownTagKeys order.
     */
    std::vector<std::string> TagSubModelPartKeys() const {
        const std::vector<std::string>& source = mTagKeys.empty() ? KnownTagKeys() : mTagKeys;
        std::vector<std::string> out;
        out.reserve(source.size());
        for (const std::string& r_key : source)
            if (std::find(mExcludedTagKeys.begin(), mExcludedTagKeys.end(), r_key) ==
                mExcludedTagKeys.end())
                out.push_back(r_key);
        return out;
    }

private:
    /** @brief Per staged block: what it became in the ModelPart. */
    struct BlockRecord {
        enum class Kind { Element, Condition, Ragged } mKind = Kind::Ragged;
        IndexType mFirstId = 0;  // first entity Id (Element/Condition kinds)
        std::size_t mCount = 0;
        /// The block's Kratos entity name, "" to derive it from the cell type.
        std::string mEntityName;
    };

    void ResetModelPartOnly() {
        // Ingestion after materialization restarts the ModelPart view.
        if (mMaterialized) {
            mpRoot.reset();
            mRecords.clear();
            mMaterialized = false;
            mStale = false;
        }
    }

    /** @brief Drop a built ModelPart so the next access rebuilds it. */
    void InvalidateMaterialization() {
        if (mMaterialized && !mStale) {
            mMaterialized = false;
            mpRoot.reset();
        }
    }

    /** @brief Shared body of the two `GetModelPart()` overloads. */
    ModelPart& EnsureMaterialized() const {
        EnsureStage();  // a pending user mutation must be folded in first
        if (!mMaterialized)
            Materialize();
        return *mpRoot;
    }

    const NativeMesh& Stage() const {
        EnsureStage();
        return mStage;
    }
    void EnsureStage() const {
        if (mStale) {
            RebuildStageFromModelPart();
            mStale = false;
        }
    }

    /** @brief Topological dimension of a staged block (Custom via the name table). */
    static int BlockDimension(const NativeCellBlock& rBlock) {
        const int dim = cell_type_dimension(rBlock.mType);
        if (dim >= 0)
            return dim;
        auto it = topological_dimension().find(rBlock.mTypeName);
        if (it != topological_dimension().end())
            return it->second;
        return 3;
    }

    void Materialize() const {
        mpRoot = std::make_unique<ModelPart>("Main");
        mRecords.clear();
        ModelPart& r_mp = *mpRoot;

        // Nodes: Id = index + 1, z zero-padded for 2-D points.
        const NDArray& points = mStage.Points();
        const std::size_t npts = mStage.NumPoints();
        const std::size_t dim = mStage.PointDim();
        const double* p = npts ? points.As<double>() : nullptr;
        for (std::size_t i = 0; i < npts; ++i)
            r_mp.CreateNewNode(i + 1, p[i * dim], dim > 1 ? p[i * dim + 1] : 0.0,
                               dim > 2 ? p[i * dim + 2] : 0.0);

        // Elements/Conditions split: block dim == mesh max dim -> Element.
        int mesh_dim = 0;
        for (const auto& r_b : mStage.Blocks())
            mesh_dim = std::max(mesh_dim, BlockDimension(r_b));

        // Per-block property ids, when the mesh carries a usable tag array.
        // Empty means "no tags", which is the pre-9.1.0 behaviour of 0 for all.
        std::vector<std::vector<std::int64_t>> prop_ids;
        if (mStage.HasCellData(kKratosPropertyTagKey) &&
            mStage.CellDataNumBlocks(kKratosPropertyTagKey) == mStage.NumCellBlocks()) {
            bool usable = true;
            for (std::size_t b = 0; b < mStage.NumCellBlocks() && usable; ++b) {
                const NDArray& r_tags = mStage.CellData(kKratosPropertyTagKey, b);
                if (detail::is_float_dtype(r_tags.Dtype()) ||
                    r_tags.Size() != mStage.Blocks()[b].NumCells()) {
                    usable = false;  // a float or mis-sized tag is not an id
                    break;
                }
                std::vector<std::int64_t> col(r_tags.Size());
                for (std::size_t i = 0; i < r_tags.Size(); ++i)
                    col[i] = detail::read_int(r_tags, i);
                prop_ids.push_back(std::move(col));
            }
            if (!usable)
                prop_ids.clear();
        }
        for (const std::vector<std::int64_t>& r_col : prop_ids)
            for (std::int64_t id : r_col)
                if (id >= 0)
                    r_mp.CreateNewProperties(static_cast<IndexType>(id));

        // Real material values, where the mesh carries them. The loop above
        // only creates the ids the entity rows reference, leaving every block
        // empty; this fills them in. CreateNewProperties is find-or-create, so
        // an id both passes see is filled in place rather than duplicated, and
        // the creation order above is unchanged. This is what makes
        // kratos_bridge's "apply property" to_model_part overload transfer
        // anything at all -- it iterates rSource.Properties(), which held
        // id-only sets before v9.2.0.
        for (const PropertySet& r_set : mStage.PropertySets().All())
            r_mp.CreateNewProperties(static_cast<IndexType>(r_set.mId)).mValues = r_set.mValues;

        IndexType next_elem = 1, next_cond = 1;
        std::size_t n_elem_rows = 0, n_cond_rows = 0;
        std::size_t block_index = 0;
        for (const auto& r_b : mStage.Blocks()) {
            BlockRecord rec;
            rec.mCount = r_b.NumCells();
            if (block_index < mBlockEntityNames.size())
                rec.mEntityName = mBlockEntityNames[block_index];
            if (r_b.IsRagged()) {
                rec.mKind = BlockRecord::Kind::Ragged;  // pass-through, no entities
            } else {
                const bool is_elem = BlockDimension(r_b) == mesh_dim;
                rec.mKind = is_elem ? BlockRecord::Kind::Element : BlockRecord::Kind::Condition;
                rec.mFirstId = is_elem ? next_elem : next_cond;
                const std::size_t n = r_b.NumCells();
                const std::size_t k = r_b.mConn.Ndim() >= 2 ? r_b.mConn.Shape()[1] : 0;
                const std::int64_t* conn = r_b.mConn.As<std::int64_t>();
                for (std::size_t c = 0; c < n; ++c) {
                    std::vector<IndexType> ids(k);
                    for (std::size_t j = 0; j < k; ++j)
                        ids[j] = static_cast<IndexType>(conn[c * k + j]) + 1;
                    // Properties id from the staged gmsh:physical tag when it
                    // covers every block: before this it was always 0, so a
                    // materialized ModelPart silently lost every material
                    // assignment the file carried.
                    const IndexType prop =
                        prop_ids.empty() ? 0 : static_cast<IndexType>(prop_ids[block_index][c]);
                    if (is_elem)
                        r_mp.CreateNewElement(r_b.mType, next_elem++, std::move(ids), prop,
                                              rec.mEntityName);
                    else
                        r_mp.CreateNewCondition(r_b.mType, next_cond++, std::move(ids), prop,
                                                rec.mEntityName);
                }
                (is_elem ? n_elem_rows : n_cond_rows) += n;
            }
            mRecords.push_back(rec);
            ++block_index;
        }

        // point_data -> nodal data (shared row order: node index).
        for (const auto& r_name : mStage.PointDataNames())
            r_mp.SetNodalData(r_name, mStage.PointData(r_name));

        // cell_data -> elemental/conditional data, concatenated per kind in
        // entity order (== block order within each kind).
        for (const auto& r_name : mStage.CellDataNames()) {
            if (mStage.CellDataNumBlocks(r_name) != mStage.NumCellBlocks())
                continue;  // partial data cannot be aligned with entities
            if (n_elem_rows > 0) {
                NDArray col = ConcatKind(r_name, BlockRecord::Kind::Element, n_elem_rows);
                if (col.Size() > 0)
                    r_mp.SetElementalData(r_name, std::move(col));
            }
            if (n_cond_rows > 0) {
                NDArray col = ConcatKind(r_name, BlockRecord::Kind::Condition, n_cond_rows);
                if (col.Size() > 0)
                    r_mp.SetConditionalData(r_name, std::move(col));
            }
        }

        // Named regions -> SubModelParts. These run FIRST and win their names:
        // a region is an explicit, format-declared group, whereas the tag pass
        // below is an inference from integer cell_data. Both go through the
        // same `HasSubModelPart` guard, so the tag pass simply skips a name a
        // region already claimed (documented in doc/regions.md).
        BuildSubModelPartsFromRegions();

        // Integer tags -> SubModelParts (unless disabled).
        if (mTagsToSubModelParts)
            for (const auto& r_key : TagSubModelPartKeys())
                if (mStage.HasCellData(r_key) &&
                    mStage.CellDataNumBlocks(r_key) == mStage.NumCellBlocks())
                    BuildSubModelPartsFor(r_key);

        mMaterialized = true;
    }

    /**
     * @brief Materialize every Point/Cell region as a named SubModelPart.
     *
     * Cell regions index cells *globally* (block-major, `detail/cell_index.hpp`);
     * `mRecords` already knows each block's first entity Id and kind, so the
     * translation to Element/Condition Ids is a scan over the block bases.
     * Point regions become node-only SubModelParts. `Side` regions have no
     * ModelPart entity to attach to and are skipped with a warning.
     */
    /** @brief Split a region name on the '/' path separator. */
    static std::vector<std::string> SplitRegionPath(const std::string& rName) {
        std::vector<std::string> out;
        std::size_t start = 0;
        while (true) {
            const std::size_t slash = rName.find('/', start);
            out.push_back(rName.substr(start, slash - start));
            if (slash == std::string::npos)
                break;
            start = slash + 1;
        }
        return out;
    }

    /** @brief Whether the SubModelPart chain @p rName names already exists. */
    bool RegionSubModelPartExists(const std::string& rName) const {
        const ModelPart* p_at = mpRoot.get();
        for (const std::string& r_seg : SplitRegionPath(rName)) {
            if (!p_at->HasSubModelPart(r_seg))
                return false;
            p_at = &p_at->GetSubModelPart(r_seg);
        }
        return true;
    }

    /**
     * @brief Walk/create the SubModelPart chain a region name denotes.
     *
     * `'/'` is the nesting separator, matching what the MDPA reader already
     * emits for a nested `Begin SubModelPart` and what `RestoreRegions` reads
     * back, so nesting survives a full region -> ModelPart -> region trip.
     * `'.'` cannot appear in a segment (ModelPart::FullName joins with it), and
     * an empty segment names nothing; either is a warning and a skip rather than
     * an exception escaping a lazy `GetModelPart()`, and the region itself stays
     * on the mesh untouched.
     *
     * @return the leaf part, or null when the name cannot be used.
     */
    ModelPart* ResolveRegionSubModelPart(const std::string& rName) const {
        const std::vector<std::string> segments = SplitRegionPath(rName);
        for (const std::string& r_seg : segments) {
            if (!r_seg.empty() && r_seg.find('.') == std::string::npos)
                continue;
            log::warn(
                "kratos backend: region '{}' has a segment that cannot name a SubModelPart "
                "(empty, or containing '.', which ModelPart::FullName reserves); it is not "
                "materialized (the region itself is kept)",
                rName);
            return nullptr;
        }
        ModelPart* p_at = mpRoot.get();
        for (const std::string& r_seg : segments)
            p_at = p_at->HasSubModelPart(r_seg) ? &p_at->GetSubModelPart(r_seg)
                                                : &p_at->CreateSubModelPart(r_seg);
        return p_at;
    }

    void BuildSubModelPartsFromRegions() const {
        const std::size_t nregions = mStage.NumRegions();
        if (nregions == 0)
            return;

        // Global cell index -> (block, row), from the staged block sizes. Built
        // once here rather than via detail::cell_index.hpp, which would need
        // mesh.hpp and so create a circular include from a backend header.
        std::vector<std::int64_t> bases;
        bases.reserve(mStage.NumCellBlocks() + 1);
        bases.push_back(0);
        for (const auto& r_b : mStage.Blocks())
            bases.push_back(bases.back() + static_cast<std::int64_t>(r_b.NumCells()));

        for (std::size_t i = 0; i < nregions; ++i) {
            const meshioplusplus::Region& r_region = mStage.Region(i);
            if (r_region.mKind == RegionKind::Side) {
                log::warn(
                    "kratos backend: region '{}' is a side set; SubModelParts have no "
                    "facet entities, so it is not materialized (the region itself is kept)",
                    r_region.mName);
                continue;
            }
            if (RegionSubModelPartExists(r_region.mName))
                continue;

            const std::int64_t* entries = r_region.Entries();
            const std::size_t n = r_region.NumEntries();
            ModelPart* p_smp = ResolveRegionSubModelPart(r_region.mName);
            if (!p_smp)
                continue;  // unusable name; warned about, and the region is kept
            ModelPart& r_smp = *p_smp;

            if (r_region.mKind == RegionKind::Point) {
                std::vector<IndexType> node_ids;
                node_ids.reserve(n);
                const auto npts = static_cast<std::int64_t>(mStage.NumPoints());
                for (std::size_t k = 0; k < n; ++k)
                    if (entries[k] >= 0 && entries[k] < npts)
                        node_ids.push_back(static_cast<IndexType>(entries[k]) + 1);
                r_smp.AddNodes(node_ids);
                continue;
            }

            std::vector<IndexType> elems, conds;
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t g = entries[k];
                if (g < 0 || bases.empty() || g >= bases.back())
                    continue;
                std::size_t b = 0;
                while (b + 1 < bases.size() && g >= bases[b + 1])
                    ++b;
                const BlockRecord& rec = mRecords[b];
                if (rec.mKind == BlockRecord::Kind::Ragged)
                    continue;  // ragged blocks have no ModelPart entities
                const IndexType id = rec.mFirstId + static_cast<IndexType>(g - bases[b]);
                (rec.mKind == BlockRecord::Kind::Element ? elems : conds).push_back(id);
            }
            r_smp.AddElements(elems);
            r_smp.AddConditions(conds);
            // Kratos convention: a sub model part contains its entities' nodes.
            std::vector<IndexType> node_ids;
            detail::IdList seen;
            for (IndexType eid : elems)
                for (IndexType nid : mpRoot->GetElement(eid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            for (IndexType cid : conds)
                for (IndexType nid : mpRoot->GetCondition(cid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            r_smp.AddNodes(node_ids);
        }
    }

    /** @brief Concatenate one cell-data name's arrays over blocks of one kind. */
    NDArray ConcatKind(const std::string& rName, BlockRecord::Kind kind,
                       std::size_t totalRows) const {
        // Determine dtype/trailing shape from the first contributing block;
        // bail out (empty result) if the blocks disagree on dtype.
        const NDArray* p_first = nullptr;
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind != kind || mRecords[b].mCount == 0)
                continue;
            const NDArray& blk = mStage.CellData(rName, b);
            if (!p_first)
                p_first = &blk;
            else if (blk.Dtype() != p_first->Dtype())
                return NDArray{};
        }
        if (!p_first || totalRows == 0)
            return NDArray{};
        std::vector<std::size_t> shape = p_first->Shape();
        if (shape.empty())
            shape = {0};
        shape[0] = totalRows;
        NDArray out = NDArray::Uninit(p_first->Dtype(), shape);
        std::size_t off = 0;
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind != kind)
                continue;
            const NDArray& blk = mStage.CellData(rName, b);
            std::memcpy(out.Data() + off, blk.Data(), blk.Nbytes());
            off += blk.Nbytes();
        }
        return out;
    }

    void BuildSubModelPartsFor(const std::string& rKey) const {
        // Only integer-kind scalar-per-cell arrays qualify as tags.
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind == BlockRecord::Kind::Ragged)
                continue;
            const NDArray& a = mStage.CellData(rKey, b);
            if (detail::is_float_dtype(a.Dtype()) || a.Size() != mRecords[b].mCount)
                return;
        }
        std::string prefix = rKey;
        for (char& r_c : prefix)
            if (r_c == ':')
                r_c = '_';

        // tag value -> (element ids, condition ids)
        struct Members {
            std::vector<IndexType> mElems, mConds;
        };
        std::unordered_map<std::int64_t, Members> groups;
        std::vector<std::int64_t> order;  // first-seen order for determinism
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            const BlockRecord& rec = mRecords[b];
            if (rec.mKind == BlockRecord::Kind::Ragged)
                continue;
            const NDArray& a = mStage.CellData(rKey, b);
            for (std::size_t c = 0; c < rec.mCount; ++c) {
                const std::int64_t tag = detail::read_int(a, c);
                auto [it, inserted] = groups.try_emplace(tag);
                if (inserted)
                    order.push_back(tag);
                if (rec.mKind == BlockRecord::Kind::Element)
                    it->second.mElems.push_back(rec.mFirstId + c);
                else
                    it->second.mConds.push_back(rec.mFirstId + c);
            }
        }
        for (const std::int64_t tag : order) {
            const std::string name = prefix + "_" + std::to_string(tag);
            if (mpRoot->HasSubModelPart(name))
                continue;  // an earlier tag key already claimed the name
            ModelPart& r_smp = mpRoot->CreateSubModelPart(name);
            const Members& r_m = groups.at(tag);
            r_smp.AddElements(r_m.mElems);
            r_smp.AddConditions(r_m.mConds);
            // Kratos convention: a sub model part contains its entities' nodes.
            std::vector<IndexType> node_ids;
            detail::IdList seen;
            for (IndexType eid : r_m.mElems)
                for (IndexType nid : mpRoot->GetElement(eid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            for (IndexType cid : r_m.mConds)
                for (IndexType nid : mpRoot->GetCondition(cid).NodeIds())
                    if (seen.Add(nid))
                        node_ids.push_back(nid);
            r_smp.AddNodes(node_ids);
        }
    }

    /**
     * @brief The name an entity would be written under: its own, or the derived
     * default for its cell type.
     */
    static std::string EffectiveEntityName(const GeometricalEntity& rEntity,
                                           BlockRecord::Kind kind) {
        if (rEntity.HasName())
            return rEntity.Name();
        return kind == BlockRecord::Kind::Element ? kratos_element_name(rEntity.Type())
                                                  : kratos_condition_name(rEntity.Type());
    }

    /** @brief One block of the rebuilt mesh, before it is added to the stage. */
    struct PendingBlock {
        BlockRecord mRecord;
        std::string mType;                                     // entity blocks
        NDArray mConn;                                         // entity blocks
        std::size_t mOldBlock = static_cast<std::size_t>(-1);  // ragged: source block
    };

    void RebuildStageFromModelPart() const {
        const ModelPart& r_mp = *mpRoot;
        NativeMesh fresh;

        // Points from nodes in container order; node Id -> 0-based index.
        const std::size_t n = r_mp.Nodes().Size();
        NDArray pts = NDArray::Uninit(DType::Float64, {n, 3});
        double* p = pts.As<double>();
        std::size_t i = 0;
        for (const Node& r_node : r_mp.Nodes()) {
            p[i * 3 + 0] = r_node.X();
            p[i * 3 + 1] = r_node.Y();
            p[i * 3 + 2] = r_node.Z();
            ++i;
        }
        fresh.AssignPoints(std::move(pts));

        // Consecutive same-type runs -> blocks (Elements first, then
        // Conditions), matching the materialization convention.
        std::vector<PendingBlock> entity_blocks;
        CollectEntityBlocks(r_mp, r_mp.Elements(), BlockRecord::Kind::Element, entity_blocks);
        CollectEntityBlocks(r_mp, r_mp.Conditions(), BlockRecord::Kind::Condition, entity_blocks);

        // Interleave the preserved ragged blocks back in. They never became
        // entities, so the user's ModelPart edits cannot have touched them and
        // carrying them through verbatim is exact. Walking the OLD records puts
        // each one back where it was, which keeps the block layout (and so the
        // per-block cell_data alignment) stable for the overwhelmingly common
        // case of a mesh whose entity structure was not restructured.
        std::vector<PendingBlock> plan;
        plan.reserve(entity_blocks.size() + mRecords.size());
        std::size_t next_entity = 0;
        for (std::size_t b = 0; b < mRecords.size(); ++b) {
            if (mRecords[b].mKind == BlockRecord::Kind::Ragged) {
                PendingBlock pb;
                pb.mRecord = mRecords[b];
                pb.mOldBlock = b;
                plan.push_back(std::move(pb));
            } else if (next_entity < entity_blocks.size()) {
                plan.push_back(std::move(entity_blocks[next_entity++]));
            }
        }
        for (; next_entity < entity_blocks.size(); ++next_entity)
            plan.push_back(std::move(entity_blocks[next_entity]));

        for (const PendingBlock& r_pb : plan) {
            if (r_pb.mOldBlock != static_cast<std::size_t>(-1))
                CopyRaggedBlock(mStage.Cells(r_pb.mOldBlock), fresh);
            else
                fresh.AddCellBlock(r_pb.mType, r_pb.mConn);
        }

        // Nodal data -> point_data; elemental/conditional -> per-block slices.
        for (const auto& r_name : r_mp.NodalDataNames())
            fresh.AddPointData(r_name, r_mp.GetNodalData(r_name));
        RestoreCellData(r_mp.ElementalDataNames(), BlockRecord::Kind::Element, r_mp, plan, fresh);
        RestoreCellData(r_mp.ConditionalDataNames(), BlockRecord::Kind::Condition, r_mp, plan,
                        fresh);
        for (const auto& r_name : mStage.FieldDataNames())
            fresh.AddFieldData(r_name, mStage.FieldData(r_name));

        RestoreRegions(r_mp, plan, fresh);
        // The ModelPart is the authority after a direct mutation, so read the
        // properties back from it rather than carrying the old stage's forward:
        // a set the user deleted must not be resurrected. They are keyed by id,
        // not by entity index, so unlike regions and cell_data slices there is
        // nothing to remap.
        for (const PropertySet& r_set : r_mp.Properties())
            fresh.AddPropertySet(r_set);

        mRecords.clear();
        for (PendingBlock& r_pb : plan)
            mRecords.push_back(r_pb.mRecord);
        // Re-align the per-block names with the rebuilt block layout; a ragged
        // block never had an entity name, so it keeps an empty one.
        mBlockEntityNames.assign(plan.size(), std::string());
        for (std::size_t b = 0; b < plan.size(); ++b)
            mBlockEntityNames[b] = plan[b].mRecord.mEntityName;
        mStage = std::move(fresh);
    }

    /** @brief Re-add a ragged (polygon/polyhedron) block verbatim. */
    static void CopyRaggedBlock(const NativeMesh::CellView& rView, NativeMesh& rOut) {
        const std::string& r_type = rView.Type();
        const std::size_t nc = rView.NumCells();
        if (rView.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(nc);
            for (std::size_t c = 0; c < nc; ++c) {
                const std::size_t nf = rView.NumFaces(c);
                cells[c].resize(nf);
                for (std::size_t f = 0; f < nf; ++f) {
                    const auto [p_face, len] = rView.Face(c, f);
                    cells[c][f].assign(p_face, p_face + len);
                }
            }
            rOut.AddPolyhedronBlock(r_type, std::move(cells));
        } else {
            std::vector<std::vector<std::int64_t>> rows(nc);
            for (std::size_t c = 0; c < nc; ++c) {
                const std::int64_t* p_row = rView.Row(c);
                rows[c].assign(p_row, p_row + rView.RowSize(c));
            }
            rOut.AddPolygonBlock(r_type, std::move(rows));
        }
    }

    template <class TContainer>
    void CollectEntityBlocks(const ModelPart& rMp, const TContainer& rEntities,
                             BlockRecord::Kind kind, std::vector<PendingBlock>& rOut) const {
        std::vector<const GeometricalEntity*> run;
        auto flush = [&]() {
            if (run.empty())
                return;
            const std::size_t nc = run.size();
            const std::size_t k = run.front()->NumberOfNodes();
            NDArray conn = NDArray::Uninit(DType::Int64, {nc, k});
            std::int64_t* c = conn.As<std::int64_t>();
            for (std::size_t r = 0; r < nc; ++r)
                for (std::size_t j = 0; j < k; ++j)
                    c[r * k + j] =
                        static_cast<std::int64_t>(rMp.Nodes().IndexOf(run[r]->NodeIds()[j]));
            PendingBlock pb;
            pb.mRecord.mKind = kind;
            pb.mRecord.mFirstId = run.front()->Id();
            pb.mRecord.mCount = nc;
            pb.mRecord.mEntityName = run.front()->Name();
            pb.mType = cell_type_name(run.front()->Type());
            pb.mConn = std::move(conn);
            rOut.push_back(std::move(pb));
            run.clear();
        };
        for (const auto& r_e : rEntities) {
            // The Kratos entity NAME is part of the run key, not just the cell
            // type: two adjacent SmallDisplacementElement3D4N and
            // TotalLagrangianElement3D4N runs are both `tetra`, and merging them
            // would leave one name to write both back under.
            //
            // The comparison is on the EFFECTIVE name -- an entity created from
            // a bare CellType has no name and means "derive it" -- so adding an
            // explicitly-named `Element3D4N` next to entities materialized
            // without names still groups into one block, as it always did.
            if (!run.empty() &&
                (run.front()->Type() != r_e.Type() ||
                 run.front()->NumberOfNodes() != r_e.NumberOfNodes() ||
                 EffectiveEntityName(*run.front(), kind) != EffectiveEntityName(r_e, kind)))
                flush();
            run.push_back(&r_e);
        }
        flush();
    }

    void RestoreCellData(const std::vector<std::string>& rNames, BlockRecord::Kind kind,
                         const ModelPart& rMp, const std::vector<PendingBlock>& rPlan,
                         NativeMesh& rOut) const {
        for (const auto& r_name : rNames) {
            const NDArray& col = kind == BlockRecord::Kind::Element
                                     ? rMp.GetElementalData(r_name)
                                     : rMp.GetConditionalData(r_name);
            const std::size_t ncols = col.Ndim() >= 2 ? col.Shape()[1] : 1;
            const std::size_t isz = dtype_size(col.Dtype());
            // A ragged block's slice never reached the ModelPart, so it can only
            // come from the old stage -- and only if that array covered every
            // block there (AppendCellData demands one slice per block, so a
            // partially-covered array cannot be reassembled at all).
            const bool old_has = mStage.HasCellData(r_name) &&
                                 mStage.CellDataNumBlocks(r_name) == mStage.NumCellBlocks();
            std::size_t row = 0;
            for (const PendingBlock& r_pb : rPlan) {
                if (r_pb.mOldBlock != static_cast<std::size_t>(-1)) {
                    if (old_has)
                        rOut.AppendCellData(r_name, mStage.CellData(r_name, r_pb.mOldBlock));
                    continue;
                }
                if (r_pb.mRecord.mKind != kind)
                    continue;
                std::vector<std::size_t> shape = col.Shape();
                if (shape.empty())
                    shape = {0};
                shape[0] = r_pb.mRecord.mCount;
                NDArray slice = NDArray::Uninit(col.Dtype(), shape);
                std::memcpy(slice.Data(), col.Data() + row * ncols * isz,
                            r_pb.mRecord.mCount * ncols * isz);
                row += r_pb.mRecord.mCount;
                rOut.AppendCellData(r_name, std::move(slice));
            }
        }
    }

    /**
     * @brief Read the SubModelPart structure back as named regions.
     *
     * A SubModelPart carries a name and its member entities but nothing else, so
     * the `mDim`/`mTag` of a reconstructed region are taken from a staged region
     * of the same name and kind when one exists, and are `-1` otherwise. A
     * staged region whose name no longer has a SubModelPart rides through
     * unchanged -- that is what keeps `Side` regions (which `Materialize()`
     * deliberately never turns into SubModelParts) alive across a rebuild.
     */
    void RestoreRegions(const ModelPart& rMp, const std::vector<PendingBlock>& rPlan,
                        NativeMesh& rOut) const {
        // Entity Id -> global block-major cell index, over the rebuilt layout.
        std::unordered_map<IndexType, std::int64_t> elem_to_cell, cond_to_cell;
        std::int64_t global = 0;
        for (const PendingBlock& r_pb : rPlan) {
            const bool ragged = r_pb.mOldBlock != static_cast<std::size_t>(-1);
            if (!ragged) {
                auto& r_map =
                    r_pb.mRecord.mKind == BlockRecord::Kind::Element ? elem_to_cell : cond_to_cell;
                for (std::size_t c = 0; c < r_pb.mRecord.mCount; ++c)
                    r_map[static_cast<IndexType>(r_pb.mRecord.mFirstId + c)] =
                        global + static_cast<std::int64_t>(c);
            }
            global += static_cast<std::int64_t>(r_pb.mRecord.mCount);
        }

        auto staged_meta = [&](const std::string& rName, RegionKind kind, int& rDim, int& rTag) {
            for (const auto& r_region : mStage.Regions().All())
                if (r_region.mName == rName && r_region.mKind == kind) {
                    rDim = r_region.mDim;
                    rTag = r_region.mTag;
                    return;
                }
            rDim = -1;
            rTag = -1;
        };

        std::vector<std::string> rebuilt;
        // Recursive, with the '/' path convention: a nested SubModelPart comes
        // back as one region named `parent/child`, which is exactly what the
        // MDPA reader already produces and what BuildSubModelPartsFromRegions
        // splits again on the way in. Before this the walk was top-level only,
        // so every nested part was silently dropped by a rebuild.
        RestoreRegionsFrom(rMp, rMp, std::string(), rPlan, elem_to_cell, cond_to_cell, staged_meta,
                           rebuilt, rOut);

        for (const auto& r_region : mStage.Regions().All())
            if (std::find(rebuilt.begin(), rebuilt.end(), r_region.mName) == rebuilt.end())
                rOut.AddRegion(r_region);
    }

    /** @brief One level of `RestoreRegions`' walk; recurses into nested parts. */
    template <class TStagedMeta>
    void RestoreRegionsFrom(const ModelPart& rMp, const ModelPart& rParent,
                            const std::string& rPrefix, const std::vector<PendingBlock>& rPlan,
                            const std::unordered_map<IndexType, std::int64_t>& rElemToCell,
                            const std::unordered_map<IndexType, std::int64_t>& rCondToCell,
                            const TStagedMeta& rStagedMeta, std::vector<std::string>& rRebuilt,
                            NativeMesh& rOut) const {
        (void)rPlan;
        for (const auto& r_name : rParent.SubModelPartNames()) {
            const ModelPart& r_smp = rParent.GetSubModelPart(r_name);
            const std::string full = rPrefix.empty() ? r_name : rPrefix + "/" + r_name;

            std::vector<std::int64_t> cells;
            for (IndexType id : r_smp.ElementIds()) {
                const auto it = rElemToCell.find(id);
                if (it != rElemToCell.end())
                    cells.push_back(it->second);
            }
            for (IndexType id : r_smp.ConditionIds()) {
                const auto it = rCondToCell.find(id);
                if (it != rCondToCell.end())
                    cells.push_back(it->second);
            }
            if (!cells.empty()) {
                int dim = -1, tag = -1;
                rStagedMeta(full, RegionKind::Cell, dim, tag);
                rOut.AddRegion(MakeRegion(full, RegionKind::Cell, dim, tag, cells));
                rRebuilt.push_back(full);
            }

            // A node-only SubModelPart is a Point region; one that also has
            // entities carries their nodes implicitly, so emitting a Point
            // region there too would invent a group the mesh never had.
            if (cells.empty()) {
                std::vector<std::int64_t> pts;
                for (IndexType id : r_smp.NodeIds())
                    pts.push_back(static_cast<std::int64_t>(rMp.Nodes().IndexOf(id)));
                if (!pts.empty()) {
                    int dim = -1, tag = -1;
                    rStagedMeta(full, RegionKind::Point, dim, tag);
                    rOut.AddRegion(MakeRegion(full, RegionKind::Point, dim, tag, pts));
                    rRebuilt.push_back(full);
                }
            }

            RestoreRegionsFrom(rMp, r_smp, full, rPlan, rElemToCell, rCondToCell, rStagedMeta,
                               rRebuilt, rOut);
        }
    }

    static meshioplusplus::Region MakeRegion(const std::string& rName, RegionKind kind, int dim,
                                             int tag, const std::vector<std::int64_t>& rEntries) {
        NDArray a = NDArray::Uninit(DType::Int64, {rEntries.size()});
        std::memcpy(a.Data(), rEntries.data(), rEntries.size() * sizeof(std::int64_t));
        meshioplusplus::Region region;
        region.mName = rName;
        region.mKind = kind;
        region.mDim = dim;
        region.mTag = tag;
        region.mEntries = std::move(a);
        return region;
    }

    mutable NativeMesh mStage;  // staging + writer-serving storage
    // mutable: the ModelPart is a lazily-built cache of the staged mesh, so
    // the const GetModelPart() overload can materialize it on first use.
    mutable std::unique_ptr<ModelPart> mpRoot;
    mutable std::vector<BlockRecord> mRecords;
    mutable bool mMaterialized = false;
    mutable bool mStale = false;
    /// Per staged block, the Kratos entity name to create its entities under.
    /// `mutable` for the same reason `mRecords` is: the const `EnsureStage()`
    /// rebuild re-derives it from the ModelPart.
    mutable std::vector<std::string> mBlockEntityNames;
    bool mTagsToSubModelParts = true;
    /// Empty = every KnownTagKeys() entry; see SetTagSubModelPartKeys.
    std::vector<std::string> mTagKeys;
    std::vector<std::string> mExcludedTagKeys;
};

}  // namespace meshioplusplus
