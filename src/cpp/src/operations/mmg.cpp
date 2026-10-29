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
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_MMG
#include <mmg/libmmg.h>
#endif

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/lsq_fit.hpp"
#include "../detail/mmg_guard.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

#ifdef MESHIOPLUSPLUS_HAS_MMG
// The width CMake read from mmgcmakedefines.h must be the one the headers use.
#ifdef MESHIOPLUSPLUS_MMG_INT64
static_assert(sizeof(MMG5_int) == sizeof(std::int64_t), "mmg: MMG5_int is not 64-bit");
#else
static_assert(sizeof(MMG5_int) == sizeof(std::int32_t), "mmg: MMG5_int is not 32-bit");
#endif
#endif

bool mmg_available() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_MMG
    return true;
#else
    return false;
#endif
}

bool mmg_has_elas() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_MMG_ELAS
    return true;
#else
    return false;
#endif
}

int mmg_index_bits() noexcept {
#if !defined(MESHIOPLUSPLUS_HAS_MMG)
    return 0;
#elif defined(MESHIOPLUSPLUS_MMG_INT64)
    return 64;
#else
    return 32;
#endif
}

namespace {

[[noreturn]] void mmg_absent(const char* pOp) {
    throw Unsupported(std::string(pOp) +
                      ": the MMG backend is not built in: rebuild with "
                      "-DMESHIOPLUSPLUS_WITH_MMG=ON and a shared MMG install "
                      "(build/build-mmg.sh; see doc/mmg.md)");
}

}  // namespace

// --- metric_from_hessian: no MMG needed ---------------------------------------------------

Mesh metric_from_hessian(const Mesh& rMesh, const MetricFromHessianOptions& rOptions) {
    const char* op = "metric_from_hessian";
    if (!rMesh.HasPointData(rOptions.mHessian))
        throw std::invalid_argument(std::string(op) + ": no point_data '" + rOptions.mHessian +
                                    "'");
    if (!(rOptions.mError > 0.0) || !(rOptions.mHmin > 0.0) || !(rOptions.mHmax >= rOptions.mHmin))
        throw std::invalid_argument(std::string(op) +
                                    ": error and hmin must be positive and hmax >= hmin");
    if (rOptions.mOutput.empty())
        throw std::invalid_argument(std::string(op) + ": the output name is empty");
    const NDArray& h = rMesh.PointData(rOptions.mHessian);
    const std::size_t n = rMesh.NumPoints();
    const std::size_t width = n ? h.Size() / n : 0;
    if (h.Size() != n * width || (width != 9 && width != 4))
        throw std::invalid_argument(std::string(op) + ": '" + rOptions.mHessian +
                                    "' needs 9 components (4 in 2-D) per point, a row-major "
                                    "Hessian such as hessian() writes");
    const std::size_t d = width == 9 ? 3 : 2;
    const double c = 0.5 * std::pow(static_cast<double>(d) / static_cast<double>(d + 1), 2.0);
    const double lo = 1.0 / (rOptions.mHmax * rOptions.mHmax);
    const double hi = 1.0 / (rOptions.mHmin * rOptions.mHmin);
    const std::size_t out_width = d == 3 ? 6 : 3;
    NDArray metric(DType::Float64, {n, out_width});
    double* m = metric.As<double>();
    const detail::DoubleView hv(h);
    for (std::size_t i = 0; i < n; ++i) {
        // The symmetrised Hessian as a 3x3 (a 2-D one padded with an inert axis).
        double a[9] = {0, 0, 0, 0, 0, 0, 0, 0, 1};
        for (std::size_t r = 0; r < d; ++r)
            for (std::size_t q = 0; q < d; ++q)
                a[r * 3 + q] = 0.5 * (hv[i * width + r * d + q] + hv[i * width + q * d + r]);
        double w[3], v[9];
        detail::lsq_sym3_eigen(a, w, v);
        double full[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (std::size_t k = 0; k < 3; ++k) {
            // The padding axis has eigenvector e3 and no say in a 2-D metric.
            if (d == 2 && std::abs(v[k * 3 + 2]) > 0.5)
                continue;
            double lambda = c * std::abs(w[k]) / rOptions.mError;
            if (!std::isfinite(lambda))
                lambda = hi;
            lambda = std::clamp(lambda, lo, hi);
            for (std::size_t r = 0; r < 3; ++r)
                for (std::size_t q = 0; q < 3; ++q)
                    full[r * 3 + q] += lambda * v[k * 3 + r] * v[k * 3 + q];
        }
        double* row = m + i * out_width;
        // The lower triangle, row by row: the order of a Medit .sol.
        row[0] = full[0];
        row[1] = full[3];
        row[2] = full[4];
        if (d == 3) {
            row[3] = full[6];
            row[4] = full[7];
            row[5] = full[8];
        }
    }
    Mesh out = detail::clone_mesh(rMesh);
    out.AddPointData(rOptions.mOutput, std::move(metric));
    return out;
}

#ifndef MESHIOPLUSPLUS_HAS_MMG

AdaptResult adapt(const Mesh&, const AdaptOptions&) {
    mmg_absent("adapt");
}

DiscretizeLevelsetResult discretize_levelset(const Mesh&, const DiscretizeLevelsetOptions&) {
    mmg_absent("discretize_levelset");
}

MeshDomainResult mesh_domain(const Mesh&, const MeshDomainOptions&) {
    mmg_absent("mesh_domain");
}

AdaptResult move(const Mesh&, const MoveOptions&) {
    mmg_absent("move");
}

#else

namespace {

enum class MmgLib { Mmg2d, Mmgs, Mmg3d };

const char* mmg_lib_name(MmgLib Lib) {
    switch (Lib) {
        case MmgLib::Mmg2d:
            return "mmg2d";
        case MmgLib::Mmgs:
            return "mmgs";
        default:
            return "mmg3d";
    }
}

/// What an MMG reference stands for: the regions and the `medit:ref` of the
/// entities that carry it.
struct MmgSignature {
    bool mFacet = false;   ///< given to MMG as a facet (MMG reuses element refs on new facets)
    bool mIsCell = false;  ///< a cell of the input (not only a side facet)
    bool mHasMeditRef = false;
    std::int64_t mMeditRef = 0;
    std::vector<std::size_t> mCellRegions;  ///< indices into the mesh's regions
    std::vector<std::size_t> mSideRegions;
};

/// The MMG references: one per distinct signature, numbered from `Base() + 1`.
class MmgRefTable {
public:
    /// Number the references from `Base + 1`, clear of the values MMG gives
    /// elements and facets itself (a level set's 2, 3 and 10).
    void SetBase(MMG5_int Base) { mBase = Base; }
    MMG5_int Ref(const MmgSignature& rSig) {
        std::vector<std::int64_t> key = {rSig.mFacet ? 1 : 0, rSig.mIsCell ? 1 : 0,
                                         rSig.mHasMeditRef ? 1 : 0, rSig.mMeditRef, -1};
        for (const std::size_t r : rSig.mCellRegions)
            key.push_back(static_cast<std::int64_t>(r));
        key.push_back(-2);
        for (const std::size_t r : rSig.mSideRegions)
            key.push_back(static_cast<std::int64_t>(r));
        auto [it, inserted] =
            mIndex.try_emplace(key, mBase + static_cast<MMG5_int>(mSigs.size() + 1));
        if (inserted)
            mSigs.push_back(rSig);
        return it->second;
    }
    /// The signature of a reference MMG returns; null for 0 or one it made up.
    const MmgSignature* Find(MMG5_int Ref) const {
        if (Ref <= mBase || static_cast<std::size_t>(Ref - mBase) > mSigs.size())
            return nullptr;
        return &mSigs[static_cast<std::size_t>(Ref - mBase - 1)];
    }
    /// The signature of a reference on an entity of that kind (facet or
    /// element): MMG gives the facets it creates the reference of an element.
    const MmgSignature* Find(MMG5_int Ref, bool Facet) const {
        const MmgSignature* s = Find(Ref);
        return s && s->mFacet == Facet ? s : nullptr;
    }
    std::size_t Size() const { return mSigs.size(); }

private:
    MMG5_int mBase = 0;
    std::map<std::vector<std::int64_t>, MMG5_int> mIndex;
    std::vector<MmgSignature> mSigs;
};

/// What an MMG operation gives MMG: a mesh to adapt, a mesh to cut by a level
/// set (its elements keep their raw `medit:ref`, which `Set_multiMat` reads),
/// the boundary lines of a planar domain to mesh, or a mesh to move (its
/// references clear of the moving facets' 10).
enum class MmgEncodeMode { Adapt, LevelSet, Domain, Move };

/// A mesh laid out for MMG: 1-based connectivity, references and attributes.
struct MmgInput {
    MmgLib mLib = MmgLib::Mmg3d;
    int mDim = 3;               ///< coordinates per vertex given to MMG
    std::size_t mPointDim = 3;  ///< the input's
    double mPlaneZ = 0.0;       ///< the z of a planar 3-D mesh given to mmg2d
    std::vector<double> mCoords;
    std::vector<MMG5_int> mVertRefs;
    // Elements: the volume (tetra, prism) or 2-D (triangle, quad) cells.
    std::vector<MMG5_int> mTets, mTetRefs, mPrisms, mPrismRefs;
    std::vector<MMG5_int> mTris, mTriRefs, mQuads, mQuadRefs;
    std::vector<MMG5_int> mEdges, mEdgeRefs;
    // 1-based ids within their kind.
    std::vector<MMG5_int> mReqVerts, mCorners, mReqTets, mReqTris, mReqEdges, mRidges;
    MmgRefTable mRefs;
    bool mHadFacetCells = false;  ///< the input had triangle (mmg3d) or line cells
    bool mHadPointRef = false;
    bool mHadCellRef = false;
    /// The elements carry their `medit:ref` (or 0) rather than a signature.
    bool mRawElementRefs = false;
    std::size_t mNumElements = 0;
    /// Input points that are required, by index, to put point regions back.
    std::vector<std::uint8_t> mRequiredPoint;
};

bool mmg_is_int(DType T) {
    return T == DType::Int8 || T == DType::Int16 || T == DType::Int32 || T == DType::Int64 ||
           T == DType::UInt8 || T == DType::UInt16 || T == DType::UInt32 || T == DType::UInt64;
}

/// The facets of the input given to MMG: sorted nodes -> its signature and nodes.
struct MmgFacet {
    std::vector<std::int64_t> mNodes;
    MmgSignature mSig;
    bool mRequired = false;
};

void mmg_add_unique(std::vector<std::size_t>& rList, std::size_t Value) {
    if (std::find(rList.begin(), rList.end(), Value) == rList.end()) {
        rList.push_back(Value);
        std::sort(rList.begin(), rList.end());
    }
}

MmgInput mmg_encode(const Mesh& rMesh, const MmgParameters& rParams, const char* pOp,
                    const std::vector<std::string>& rConsumed, std::vector<std::string>& rWarnings,
                    MmgEncodeMode Mode = MmgEncodeMode::Adapt) {
    const std::string op(pOp);
    MmgInput in;
    if (Mode == MmgEncodeMode::LevelSet) {
        in.mRawElementRefs = true;
        in.mRefs.SetBase(100000);
    } else if (Mode == MmgEncodeMode::Move) {
        in.mRefs.SetBase(100000);  // clear of kMmgDispRef, the moving facets' 10
    }
    const std::size_t np = rMesh.NumPoints();
    in.mPointDim = rMesh.PointDim();
    if (in.mPointDim != 2 && in.mPointDim != 3)
        throw std::invalid_argument(op + ": the points must be 2-D or 3-D");

    // Which library: tetrahedra go to mmg3d; triangles to mmg2d when planar.
    bool has_tet = false, has_tri = false, has_quad = false, has_wedge = false;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        if (cb.IsRagged())
            throw std::invalid_argument(op + ": a ragged " + type +
                                        " block; convert_cells it to simplices first");
        if (type == "tetra")
            has_tet = true;
        else if (type == "triangle")
            has_tri = true;
        else if (type == "quad")
            has_quad = true;
        else if (type == "wedge")
            has_wedge = true;
        else if (type != "line" && type != "vertex")
            throw std::invalid_argument(op +
                                        ": MMG takes tetrahedra, prisms, triangles, quads "
                                        "and lines, not " +
                                        type + "; convert_cells(simplexify) it first");
    }
    const detail::DoubleView pts(rMesh.Points());
    if (Mode == MmgEncodeMode::Domain) {
        if (has_tet || has_wedge || has_tri || has_quad)
            throw std::invalid_argument(op +
                                        ": give the boundary as line cells only; a mesh with "
                                        "triangles or tetrahedra is adapt's");
        bool has_line = false;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            has_line = has_line || rMesh.Cells(b).Type() == "line";
        if (!has_line)
            throw std::invalid_argument(op + ": the boundary has no line cells");
        has_tri = true;  // mmg2d, if the points are planar (checked below)
    }
    if (has_wedge && !has_tet)
        throw std::invalid_argument(op +
                                    ": mmg3d needs tetrahedra beside prisms; "
                                    "convert_cells(simplexify) the mesh first");
    if (has_quad && !has_tri && !has_tet)
        throw std::invalid_argument(op +
                                    ": mmg2d needs triangles beside quads; "
                                    "convert_cells(simplexify) the mesh first");
    if (has_tet || has_wedge) {
        if (in.mPointDim != 3)
            throw std::invalid_argument(op + ": tetrahedra need 3-D points");
        in.mLib = MmgLib::Mmg3d;
    } else if (has_tri || has_quad) {
        bool planar = in.mPointDim == 2;
        if (!planar && np > 0) {
            double zmin = pts[2], zmax = pts[2], span = 0.0;
            for (std::size_t i = 0; i < np; ++i) {
                zmin = std::min(zmin, pts[i * 3 + 2]);
                zmax = std::max(zmax, pts[i * 3 + 2]);
                for (std::size_t c = 0; c < 2; ++c)
                    span = std::max(span, std::abs(pts[i * 3 + c] - pts[c]));
            }
            planar = zmax - zmin <= 1e-12 * std::max(span, 1.0);
            in.mPlaneZ = zmin;
        }
        in.mLib = planar ? MmgLib::Mmg2d : MmgLib::Mmgs;
        if (Mode == MmgEncodeMode::Domain && !planar)
            throw std::invalid_argument(op +
                                        ": the boundary must be planar (2-D points, or 3-D "
                                        "points on one z)");
        if (has_quad && in.mLib == MmgLib::Mmgs)
            throw std::invalid_argument(op +
                                        ": mmgs takes no quads; convert_cells(simplexify) "
                                        "the surface first");
    } else {
        throw std::invalid_argument(op + ": the mesh has no triangle or tetrahedron");
    }
    in.mDim = in.mLib == MmgLib::Mmg2d ? 2 : 3;

    // The index width MMG was built with (roadmap 3.2.5.2).
    const bool int64 = mmg_index_bits() == 64;
    detail::mmg_check_index_range(np, int64, "points", pOp);
    std::size_t total_cells = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        total_cells += rMesh.Cells(b).NumCells();
    detail::mmg_check_index_range(total_cells, int64, "cells", pOp);

    // Points and their references.
    in.mCoords.resize(np * static_cast<std::size_t>(in.mDim));
    for (std::size_t i = 0; i < np; ++i)
        for (int c = 0; c < in.mDim; ++c)
            in.mCoords[i * in.mDim + c] = pts[i * in.mPointDim + c];
    in.mVertRefs.assign(np, 0);
    if (rMesh.HasPointData("medit:ref") && mmg_is_int(rMesh.PointData("medit:ref").Dtype())) {
        const detail::Int64View refs(rMesh.PointData("medit:ref"));
        for (std::size_t i = 0; i < np; ++i)
            in.mVertRefs[i] = static_cast<MMG5_int>(refs[i]);
        in.mHadPointRef = true;
    }
    const auto consumed = [&](const std::string& rName) {
        return std::find(rConsumed.begin(), rConsumed.end(), rName) != rConsumed.end();
    };
    for (const std::string& name : rMesh.PointDataNames())
        if (name != "medit:ref" && !consumed(name))
            rWarnings.push_back(op + ": point_data '" + name +
                                "' is dropped: the points are new (compose with interpolate)");
    for (const std::string& name : rMesh.CellDataNames())
        if (name != "medit:ref")
            rWarnings.push_back(op + ": cell_data '" + name +
                                "' is dropped: the cells are new (compose with interpolate)");

    // Required points: the frozen mask, and required regions below.
    in.mRequiredPoint.assign(np, 0);
    if (!rParams.mFrozen.empty()) {
        if (rParams.mFrozen.size() != np)
            throw std::invalid_argument(op + ": the frozen mask has " +
                                        std::to_string(rParams.mFrozen.size()) + " entries for " +
                                        std::to_string(np) + " points");
        for (std::size_t i = 0; i < np; ++i)
            in.mRequiredPoint[i] = rParams.mFrozen[i] != 0;
    }

    // Each cell's signature: its medit:ref and its cell regions.
    std::vector<std::size_t> block_start(rMesh.NumCellBlocks() + 1, 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        block_start[b + 1] = block_start[b] + rMesh.Cells(b).NumCells();
    std::vector<MmgSignature> cell_sig(total_cells);
    const bool cell_ref = rMesh.HasCellData("medit:ref") &&
                          rMesh.CellDataNumBlocks("medit:ref") == rMesh.NumCellBlocks() &&
                          rMesh.NumCellBlocks() > 0 &&
                          mmg_is_int(rMesh.CellData("medit:ref", 0).Dtype());
    in.mHadCellRef = cell_ref;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        std::vector<std::int64_t> refs(rMesh.Cells(b).NumCells(), 0);
        if (cell_ref) {
            const detail::Int64View v(rMesh.CellData("medit:ref", b));
            for (std::size_t i = 0; i < refs.size(); ++i)
                refs[i] = v[i];
        }
        for (std::size_t i = 0; i < refs.size(); ++i) {
            MmgSignature& s = cell_sig[block_start[b] + i];
            s.mIsCell = true;
            s.mHasMeditRef = cell_ref;
            s.mMeditRef = refs[i];
        }
    }
    std::vector<std::uint8_t> cell_required(total_cells, 0);
    std::vector<std::uint8_t> point_corner(np, 0);
    const auto required = [&](const std::string& rName) {
        return std::find(rParams.mRequired.begin(), rParams.mRequired.end(), rName) !=
               rParams.mRequired.end();
    };
    for (const std::string& name : rParams.mRequired)
        if (!rMesh.HasRegion(name))
            throw std::invalid_argument(op + ": no region '" + name + "' to require");
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        const std::int64_t* e = reg.Entries();
        const std::size_t ne = reg.NumEntries();
        const bool req = required(reg.mName) || reg.mName == "medit:RequiredVertices" ||
                         reg.mName == "medit:RequiredEdges" ||
                         reg.mName == "medit:RequiredTriangles";
        if (reg.mKind == RegionKind::Point) {
            for (std::size_t k = 0; k < ne; ++k) {
                if (e[k] < 0 || static_cast<std::size_t>(e[k]) >= np)
                    continue;
                if (req || reg.mName == "medit:Corners")
                    in.mRequiredPoint[static_cast<std::size_t>(e[k])] = 1;
                if (reg.mName == "medit:Corners")
                    point_corner[static_cast<std::size_t>(e[k])] = 1;
            }
        } else if (reg.mKind == RegionKind::Cell) {
            for (std::size_t k = 0; k < ne; ++k) {
                if (e[k] < 0 || static_cast<std::size_t>(e[k]) >= total_cells)
                    continue;
                mmg_add_unique(cell_sig[static_cast<std::size_t>(e[k])].mCellRegions, r);
                if (req)
                    cell_required[static_cast<std::size_t>(e[k])] = 1;
            }
        }
    }

    // The facets MMG is given: the input's own (triangles and quads in a
    // volume, lines elsewhere) merged with the sides of side regions.
    std::map<std::vector<std::int64_t>, MmgFacet> facets;
    const auto facet_key = [](std::vector<std::int64_t> Nodes) {
        std::sort(Nodes.begin(), Nodes.end());
        return Nodes;
    };
    std::size_t dropped_sides = 0;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind != RegionKind::Side)
            continue;
        const std::int64_t* e = reg.Entries();
        const bool req = required(reg.mName);
        for (std::size_t k = 0; k + 1 < reg.NumEntries() * 2; k += 2) {
            CellType type;
            std::vector<std::int64_t> nodes;
            if (!detail::facet_nodes(rMesh, e[k], e[k + 1], type, nodes) ||
                (type != CellType::Triangle && type != CellType::Line && type != CellType::Quad)) {
                ++dropped_sides;
                continue;
            }
            // A side of the elements MMG remeshes: faces in 3-D, edges in 2-D.
            const std::size_t cell = static_cast<std::size_t>(e[k]);
            std::size_t b = 0;
            while (b + 1 < block_start.size() && block_start[b + 1] <= cell)
                ++b;
            const std::string owner(rMesh.Cells(b).Type());
            const bool top = in.mLib == MmgLib::Mmg3d ? (owner == "tetra" || owner == "wedge")
                                                      : (owner == "triangle" || owner == "quad");
            if (!top) {
                ++dropped_sides;
                continue;
            }
            MmgFacet& f = facets[facet_key(nodes)];
            if (f.mNodes.empty())
                f.mNodes = nodes;
            mmg_add_unique(f.mSig.mSideRegions, r);
            f.mRequired = f.mRequired || req;
        }
    }
    if (dropped_sides)
        rWarnings.push_back(op + ": " + std::to_string(dropped_sides) +
                            " side entries are not on the elements MMG remeshes; dropped");

    // Elements and facets.
    const auto add_conn = [](std::vector<MMG5_int>& rOut, const detail::Int64View& rConn,
                             std::size_t Row, std::size_t K) {
        for (std::size_t j = 0; j < K; ++j)
            rOut.push_back(static_cast<MMG5_int>(rConn[Row * K + j] + 1));
    };
    // An element's reference: its signature's, or its raw medit:ref for a level set.
    const auto element_ref = [&](const MmgSignature& rSig) {
        return in.mRawElementRefs ? static_cast<MMG5_int>(rSig.mMeditRef) : in.mRefs.Ref(rSig);
    };
    std::size_t skipped_vertices = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const std::size_t k = cb.NodesPerCell();
        const detail::Int64View conn(cb.Conn());
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            const std::size_t g = block_start[b] + i;
            const MmgSignature& sig = cell_sig[g];
            if (type == "vertex") {
                ++skipped_vertices;
            } else if (type == "tetra") {
                add_conn(in.mTets, conn, i, k);
                in.mTetRefs.push_back(element_ref(sig));
                if (cell_required[g])
                    in.mReqTets.push_back(static_cast<MMG5_int>(in.mTetRefs.size()));
            } else if (type == "wedge") {
                add_conn(in.mPrisms, conn, i, k);
                in.mPrismRefs.push_back(element_ref(sig));
            } else if ((type == "triangle" || type == "quad") && in.mLib == MmgLib::Mmg3d) {
                // A boundary face: merged with the side facets on the same nodes.
                std::vector<std::int64_t> nodes(k);
                for (std::size_t j = 0; j < k; ++j)
                    nodes[j] = conn[i * k + j];
                MmgFacet& f = facets[facet_key(nodes)];
                if (f.mNodes.empty())
                    f.mNodes = nodes;
                f.mSig.mIsCell = true;
                f.mSig.mHasMeditRef = sig.mHasMeditRef;
                f.mSig.mMeditRef = sig.mMeditRef;
                for (const std::size_t r : sig.mCellRegions)
                    mmg_add_unique(f.mSig.mCellRegions, r);
                f.mRequired = f.mRequired || cell_required[g];
                in.mHadFacetCells = true;
            } else if (type == "triangle") {
                add_conn(in.mTris, conn, i, k);
                in.mTriRefs.push_back(element_ref(sig));
                if (cell_required[g])
                    in.mReqTris.push_back(static_cast<MMG5_int>(in.mTriRefs.size()));
            } else if (type == "quad") {
                add_conn(in.mQuads, conn, i, k);
                in.mQuadRefs.push_back(element_ref(sig));
            } else if (type == "line") {
                std::vector<std::int64_t> nodes(k);
                for (std::size_t j = 0; j < k; ++j)
                    nodes[j] = conn[i * k + j];
                MmgFacet& f = facets[facet_key(nodes)];
                if (f.mNodes.empty())
                    f.mNodes = nodes;
                f.mSig.mIsCell = true;
                f.mSig.mHasMeditRef = sig.mHasMeditRef;
                f.mSig.mMeditRef = sig.mMeditRef;
                for (const std::size_t r : sig.mCellRegions) {
                    mmg_add_unique(f.mSig.mCellRegions, r);
                    // A domain's boundary regions are also the triangles' sides.
                    if (Mode == MmgEncodeMode::Domain)
                        mmg_add_unique(f.mSig.mSideRegions, r);
                }
                f.mRequired = f.mRequired || cell_required[g];
                in.mHadFacetCells = true;
            }
        }
    }
    if (skipped_vertices)
        rWarnings.push_back(op + ": " + std::to_string(skipped_vertices) +
                            " vertex cells are dropped: MMG has no such element");
    for (const auto& [key, f] : facets) {
        MmgSignature sig = f.mSig;
        sig.mFacet = true;
        const MMG5_int ref = in.mRefs.Ref(sig);
        if (f.mNodes.size() == 2) {
            for (const std::int64_t v : f.mNodes)
                in.mEdges.push_back(static_cast<MMG5_int>(v + 1));
            in.mEdgeRefs.push_back(ref);
            const auto id = static_cast<MMG5_int>(in.mEdgeRefs.size());
            if (f.mRequired)
                in.mReqEdges.push_back(id);
            if (in.mLib != MmgLib::Mmg2d && f.mSig.mIsCell)
                in.mRidges.push_back(id);
        } else if (f.mNodes.size() == 3) {
            for (const std::int64_t v : f.mNodes)
                in.mTris.push_back(static_cast<MMG5_int>(v + 1));
            in.mTriRefs.push_back(ref);
            if (f.mRequired)
                in.mReqTris.push_back(static_cast<MMG5_int>(in.mTriRefs.size()));
        } else {
            for (const std::int64_t v : f.mNodes)
                in.mQuads.push_back(static_cast<MMG5_int>(v + 1));
            in.mQuadRefs.push_back(ref);
        }
    }
    for (std::size_t i = 0; i < np; ++i) {
        if (in.mRequiredPoint[i])
            in.mReqVerts.push_back(static_cast<MMG5_int>(i + 1));
        if (point_corner[i])
            in.mCorners.push_back(static_cast<MMG5_int>(i + 1));
    }
    in.mNumElements = in.mLib == MmgLib::Mmg3d ? in.mTetRefs.size() + in.mPrismRefs.size()
                                               : in.mTriRefs.size() + in.mQuadRefs.size();
    return in;
}

/// MMG's mesh and solutions, freed on every path.
struct MmgSession {
    MmgLib mLib;
    MMG5_pMesh mMesh = nullptr;
    MMG5_pSol mMet = nullptr;
    MMG5_pSol mLs = nullptr;
    MMG5_pSol mDisp = nullptr;

    MmgSession(MmgLib Lib, bool LevelSet, bool Displacement) : mLib(Lib) {
        int ok = 0;
        switch (Lib) {
            case MmgLib::Mmg3d:
                if (LevelSet)
                    ok = MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else if (Displacement)
                    ok = MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_ppDisp, &mDisp, MMG5_ARG_end);
                else
                    ok = MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_end);
                break;
            case MmgLib::Mmgs:
                if (LevelSet)
                    ok = MMGS_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                        &mMet, MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else
                    ok = MMGS_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                        &mMet, MMG5_ARG_end);
                break;
            case MmgLib::Mmg2d:
                if (LevelSet)
                    ok = MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else if (Displacement)
                    ok = MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_ppDisp, &mDisp, MMG5_ARG_end);
                else
                    ok = MMG2D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet,
                                         &mMet, MMG5_ARG_end);
                break;
        }
        if (!ok)
            throw std::runtime_error("mmg: could not initialise MMG's mesh");
    }
    ~MmgSession() {
        switch (mLib) {
            case MmgLib::Mmg3d:
                if (mLs)
                    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else if (mDisp)
                    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_ppDisp, &mDisp, MMG5_ARG_end);
                else
                    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_end);
                break;
            case MmgLib::Mmgs:
                if (mLs)
                    MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                  MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else
                    MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                  MMG5_ARG_end);
                break;
            case MmgLib::Mmg2d:
                if (mLs)
                    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_ppLs, &mLs, MMG5_ARG_end);
                else if (mDisp)
                    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_ppDisp, &mDisp, MMG5_ARG_end);
                else
                    MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mMesh, MMG5_ARG_ppMet, &mMet,
                                   MMG5_ARG_end);
                break;
        }
    }
    MmgSession(const MmgSession&) = delete;
    MmgSession& operator=(const MmgSession&) = delete;
};

void mmg_check(int Ok, const std::string& rWhat) {
    if (!Ok)
        throw std::runtime_error("mmg: MMG refused " + rWhat);
}

MMG5_int* mmg_ptr(std::vector<MMG5_int>& rV) {
    return rV.empty() ? nullptr : rV.data();
}

/// Hand the encoded mesh to MMG.
void mmg_load(MmgSession& rS, MmgInput& rIn) {
    const auto np = static_cast<MMG5_int>(rIn.mVertRefs.size());
    const auto ntet = static_cast<MMG5_int>(rIn.mTetRefs.size());
    const auto nprism = static_cast<MMG5_int>(rIn.mPrismRefs.size());
    const auto ntri = static_cast<MMG5_int>(rIn.mTriRefs.size());
    const auto nquad = static_cast<MMG5_int>(rIn.mQuadRefs.size());
    const auto nedge = static_cast<MMG5_int>(rIn.mEdgeRefs.size());
    MMG5_pMesh m = rS.mMesh;
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            mmg_check(MMG3D_Set_meshSize(m, np, ntet, nprism, ntri, nquad, nedge), "the mesh size");
            mmg_check(MMG3D_Set_vertices(m, rIn.mCoords.data(), mmg_ptr(rIn.mVertRefs)),
                      "the vertices");
            if (ntet)
                mmg_check(MMG3D_Set_tetrahedra(m, mmg_ptr(rIn.mTets), mmg_ptr(rIn.mTetRefs)),
                          "the tetrahedra");
            if (nprism)
                mmg_check(MMG3D_Set_prisms(m, mmg_ptr(rIn.mPrisms), mmg_ptr(rIn.mPrismRefs)),
                          "the prisms");
            if (ntri)
                mmg_check(MMG3D_Set_triangles(m, mmg_ptr(rIn.mTris), mmg_ptr(rIn.mTriRefs)),
                          "the triangles");
            if (nquad)
                mmg_check(MMG3D_Set_quadrilaterals(m, mmg_ptr(rIn.mQuads), mmg_ptr(rIn.mQuadRefs)),
                          "the quadrilaterals");
            if (nedge)
                mmg_check(MMG3D_Set_edges(m, mmg_ptr(rIn.mEdges), mmg_ptr(rIn.mEdgeRefs)),
                          "the edges");
            for (const MMG5_int k : rIn.mReqVerts)
                mmg_check(MMG3D_Set_requiredVertex(m, k), "a required vertex");
            for (const MMG5_int k : rIn.mCorners)
                mmg_check(MMG3D_Set_corner(m, k), "a corner");
            for (const MMG5_int k : rIn.mReqTets)
                mmg_check(MMG3D_Set_requiredTetrahedron(m, k), "a required tetrahedron");
            for (const MMG5_int k : rIn.mReqTris)
                mmg_check(MMG3D_Set_requiredTriangle(m, k), "a required triangle");
            for (const MMG5_int k : rIn.mReqEdges)
                mmg_check(MMG3D_Set_requiredEdge(m, k), "a required edge");
            for (const MMG5_int k : rIn.mRidges)
                mmg_check(MMG3D_Set_ridge(m, k), "a ridge");
            break;
        case MmgLib::Mmgs:
            mmg_check(MMGS_Set_meshSize(m, np, ntri, nedge), "the mesh size");
            mmg_check(MMGS_Set_vertices(m, rIn.mCoords.data(), mmg_ptr(rIn.mVertRefs)),
                      "the vertices");
            mmg_check(MMGS_Set_triangles(m, mmg_ptr(rIn.mTris), mmg_ptr(rIn.mTriRefs)),
                      "the triangles");
            if (nedge)
                mmg_check(MMGS_Set_edges(m, mmg_ptr(rIn.mEdges), mmg_ptr(rIn.mEdgeRefs)),
                          "the edges");
            for (const MMG5_int k : rIn.mReqVerts)
                mmg_check(MMGS_Set_requiredVertex(m, k), "a required vertex");
            for (const MMG5_int k : rIn.mCorners)
                mmg_check(MMGS_Set_corner(m, k), "a corner");
            for (const MMG5_int k : rIn.mReqTris)
                mmg_check(MMGS_Set_requiredTriangle(m, k), "a required triangle");
            for (const MMG5_int k : rIn.mReqEdges)
                mmg_check(MMGS_Set_requiredEdge(m, k), "a required edge");
            for (const MMG5_int k : rIn.mRidges)
                mmg_check(MMGS_Set_ridge(m, k), "a ridge");
            break;
        case MmgLib::Mmg2d:
            mmg_check(MMG2D_Set_meshSize(m, np, ntri, nquad, nedge), "the mesh size");
            mmg_check(MMG2D_Set_vertices(m, rIn.mCoords.data(), mmg_ptr(rIn.mVertRefs)),
                      "the vertices");
            if (ntri)
                mmg_check(MMG2D_Set_triangles(m, mmg_ptr(rIn.mTris), mmg_ptr(rIn.mTriRefs)),
                          "the triangles");
            if (nquad)
                mmg_check(MMG2D_Set_quadrilaterals(m, mmg_ptr(rIn.mQuads), mmg_ptr(rIn.mQuadRefs)),
                          "the quadrilaterals");
            if (nedge)
                mmg_check(MMG2D_Set_edges(m, mmg_ptr(rIn.mEdges), mmg_ptr(rIn.mEdgeRefs)),
                          "the edges");
            for (const MMG5_int k : rIn.mReqVerts)
                mmg_check(MMG2D_Set_requiredVertex(m, k), "a required vertex");
            for (const MMG5_int k : rIn.mCorners)
                mmg_check(MMG2D_Set_corner(m, k), "a corner");
            for (const MMG5_int k : rIn.mReqTris)
                mmg_check(MMG2D_Set_requiredTriangle(m, k), "a required triangle");
            for (const MMG5_int k : rIn.mReqEdges)
                mmg_check(MMG2D_Set_requiredEdge(m, k), "a required edge");
            break;
    }
}

int mmg_set_i(MmgSession& rS, MMG5_pSol pSol, int P3, int Ps, int P2, MMG5_int Value) {
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            return MMG3D_Set_iparameter(rS.mMesh, pSol, P3, Value);
        case MmgLib::Mmgs:
            return Ps < 0 ? 1 : MMGS_Set_iparameter(rS.mMesh, pSol, Ps, Value);
        default:
            return P2 < 0 ? 1 : MMG2D_Set_iparameter(rS.mMesh, pSol, P2, Value);
    }
}

int mmg_set_d(MmgSession& rS, MMG5_pSol pSol, int P3, int Ps, int P2, double Value) {
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            return MMG3D_Set_dparameter(rS.mMesh, pSol, P3, Value);
        case MmgLib::Mmgs:
            return MMGS_Set_dparameter(rS.mMesh, pSol, Ps, Value);
        default:
            return MMG2D_Set_dparameter(rS.mMesh, pSol, P2, Value);
    }
}

#define MIO_MMG_I(name, value)                                                                     \
    mmg_check(mmg_set_i(rS, rS.mMet, MMG3D_IPARAM_##name, MMGS_IPARAM_##name, MMG2D_IPARAM_##name, \
                        static_cast<MMG5_int>(value)),                                             \
              "the parameter " #name)
#define MIO_MMG_D(name, value)                                                                     \
    mmg_check(mmg_set_d(rS, rS.mMet, MMG3D_DPARAM_##name, MMGS_DPARAM_##name, MMG2D_DPARAM_##name, \
                        value),                                                                    \
              "the parameter " #name)

/// The shared parameters, and the per-region sizes through the reference table.
void mmg_parameters(MmgSession& rS, const MmgInput& rIn, const Mesh& rMesh, const MmgParameters& rP,
                    const char* pOp) {
    MIO_MMG_I(verbose, rP.mVerbose);
    if (rP.mMem > 0)
        MIO_MMG_I(mem, rP.mMem);
    if (rP.mHmin >= 0)
        MIO_MMG_D(hmin, rP.mHmin);
    if (rP.mHmax >= 0)
        MIO_MMG_D(hmax, rP.mHmax);
    if (rP.mHausd >= 0)
        MIO_MMG_D(hausd, rP.mHausd);
    if (rP.mHgrad >= 0)
        MIO_MMG_D(hgrad, rP.mHgrad);
    if (rP.mHgradreq >= 0)
        MIO_MMG_D(hgradreq, rP.mHgradreq);
    if (rP.mAngleDetection > 0) {
        MIO_MMG_I(angle, 1);
        MIO_MMG_D(angleDetection, rP.mAngleDetection);
    } else {
        MIO_MMG_I(angle, 0);
    }
    MIO_MMG_I(optim, rP.mOptim ? 1 : 0);
    MIO_MMG_I(noinsert, rP.mNoInsert ? 1 : 0);
    MIO_MMG_I(noswap, rP.mNoSwap ? 1 : 0);
    MIO_MMG_I(nomove, rP.mNoMove ? 1 : 0);
    MIO_MMG_I(nreg, rP.mNreg ? 1 : 0);
    if (rS.mLib != MmgLib::Mmgs) {
        if (rS.mLib == MmgLib::Mmg3d) {
            mmg_check(MMG3D_Set_iparameter(rS.mMesh, rS.mMet, MMG3D_IPARAM_nosurf, rP.mNoSurf),
                      "the parameter nosurf");
            mmg_check(MMG3D_Set_iparameter(rS.mMesh, rS.mMet, MMG3D_IPARAM_opnbdy, rP.mOpnbdy),
                      "the parameter opnbdy");
        } else {
            mmg_check(MMG2D_Set_iparameter(rS.mMesh, rS.mMet, MMG2D_IPARAM_nosurf, rP.mNoSurf),
                      "the parameter nosurf");
            mmg_check(MMG2D_Set_iparameter(rS.mMesh, rS.mMet, MMG2D_IPARAM_opnbdy, rP.mOpnbdy),
                      "the parameter opnbdy");
        }
    } else if (rP.mNoSurf || rP.mOpnbdy) {
        throw std::invalid_argument(std::string(pOp) + ": mmgs has no nosurf or opnbdy");
    }

    // Per-region sizes: every reference whose signature holds the region.
    struct Local {
        int mTyp;
        MMG5_int mRef;
        const MmgLocalParameter* pP;
    };
    std::vector<Local> locals;
    for (const MmgLocalParameter& lp : rP.mLocal) {
        std::vector<std::size_t> regions;
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
            if (rMesh.Region(r).mName == lp.mRegion)
                regions.push_back(r);
        if (regions.empty())
            throw std::invalid_argument(std::string(pOp) + ": no region '" + lp.mRegion +
                                        "' for a local parameter");
        const auto holds = [&](const std::vector<std::size_t>& rList) {
            return std::any_of(regions.begin(), regions.end(), [&](std::size_t R) {
                return std::find(rList.begin(), rList.end(), R) != rList.end();
            });
        };
        const auto add = [&](const std::vector<MMG5_int>& rRefs, int Typ) {
            std::vector<MMG5_int> seen;
            for (const MMG5_int ref : rRefs) {
                const MmgSignature* s = rIn.mRefs.Find(ref);
                if (!s || std::find(seen.begin(), seen.end(), ref) != seen.end())
                    continue;
                if (holds(s->mCellRegions) || holds(s->mSideRegions)) {
                    seen.push_back(ref);
                    locals.push_back({Typ, ref, &lp});
                }
            }
        };
        if (rS.mLib == MmgLib::Mmg3d) {
            add(rIn.mTetRefs, MMG5_Tetrahedron);
            add(rIn.mTriRefs, MMG5_Triangle);
        } else {
            add(rIn.mTriRefs, MMG5_Triangle);
            if (rS.mLib == MmgLib::Mmg2d)
                add(rIn.mEdgeRefs, MMG5_Edg);
        }
    }
    if (locals.empty())
        return;
    MIO_MMG_I(numberOfLocalParam, locals.size());
    const auto value = [](double V, double Global) { return V >= 0 ? V : Global; };
    // MMG takes all three for each: fill an unset one from the global value or
    // a value that leaves it free.
    const double ghmin = rP.mHmin >= 0 ? rP.mHmin : 0.0;
    const double ghmax = rP.mHmax >= 0 ? rP.mHmax : std::numeric_limits<double>::max() / 4;
    const double ghausd = rP.mHausd >= 0 ? rP.mHausd : 0.01;
    for (const Local& l : locals) {
        const double hmin = value(l.pP->mHmin, ghmin);
        const double hmax = value(l.pP->mHmax, ghmax);
        const double hausd = value(l.pP->mHausd, ghausd);
        int ok = 0;
        switch (rS.mLib) {
            case MmgLib::Mmg3d:
                ok = MMG3D_Set_localParameter(rS.mMesh, rS.mMet, l.mTyp, l.mRef, hmin, hmax, hausd);
                break;
            case MmgLib::Mmgs:
                ok = MMGS_Set_localParameter(rS.mMesh, rS.mMet, l.mTyp, l.mRef, hmin, hmax, hausd);
                break;
            case MmgLib::Mmg2d:
                ok = MMG2D_Set_localParameter(rS.mMesh, rS.mMet, l.mTyp, l.mRef, hmin, hmax, hausd);
                break;
        }
        mmg_check(ok, "a local parameter");
    }
}

#undef MIO_MMG_I
#undef MIO_MMG_D

/// Fill a vertex solution: `Width` values per point, `Type` MMG5_Scalar or MMG5_Tensor.
void mmg_vertex_sol(MmgSession& rS, MMG5_pSol pSol, int Type, std::vector<double>& rValues) {
    const auto np = static_cast<MMG5_int>(rS.mMesh->np);
    int ok = 0;
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            ok = MMG3D_Set_solSize(rS.mMesh, pSol, MMG5_Vertex, np, Type);
            if (ok)
                ok = Type == MMG5_Scalar   ? MMG3D_Set_scalarSols(pSol, rValues.data())
                     : Type == MMG5_Vector ? MMG3D_Set_vectorSols(pSol, rValues.data())
                                           : MMG3D_Set_tensorSols(pSol, rValues.data());
            break;
        case MmgLib::Mmgs:
            ok = MMGS_Set_solSize(rS.mMesh, pSol, MMG5_Vertex, np, Type);
            if (ok)
                ok = Type == MMG5_Scalar   ? MMGS_Set_scalarSols(pSol, rValues.data())
                     : Type == MMG5_Vector ? MMGS_Set_vectorSols(pSol, rValues.data())
                                           : MMGS_Set_tensorSols(pSol, rValues.data());
            break;
        case MmgLib::Mmg2d:
            ok = MMG2D_Set_solSize(rS.mMesh, pSol, MMG5_Vertex, np, Type);
            if (ok && Type == MMG5_Vector) {
                // MMG 5.8's MMG2D_Set_vectorSol(s) store vertex k (1-based) at
                // m[2k-1], one slot before where mmg2d's motion and its .sol
                // reader read it (m[2k]): write MMG's layout directly.
                for (MMG5_int k = 1; k <= np; ++k) {
                    pSol->m[2 * k] = rValues[static_cast<std::size_t>(2 * (k - 1))];
                    pSol->m[2 * k + 1] = rValues[static_cast<std::size_t>(2 * (k - 1) + 1)];
                }
            } else if (ok) {
                ok = Type == MMG5_Scalar ? MMG2D_Set_scalarSols(pSol, rValues.data())
                                         : MMG2D_Set_tensorSols(pSol, rValues.data());
            }
            break;
    }
    mmg_check(ok, "a solution");
}

/// A scalar point field, one value per point.
std::vector<double> mmg_scalar_field(const Mesh& rMesh, const std::string& rName, const char* pOp) {
    if (!rMesh.HasPointData(rName))
        throw std::invalid_argument(std::string(pOp) + ": no point_data '" + rName + "'");
    const NDArray& a = rMesh.PointData(rName);
    if (a.Size() != rMesh.NumPoints())
        throw std::invalid_argument(std::string(pOp) + ": '" + rName +
                                    "' must have one component per point");
    const detail::DoubleView v(a);
    std::vector<double> out(a.Size());
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = v[i];
    return out;
}

/// A metric for MMG's API (upper triangle, row by row) from the `.sol`-order field.
std::vector<double> mmg_metric_field(const Mesh& rMesh, const std::string& rName, int Dim,
                                     const char* pOp) {
    if (!rMesh.HasPointData(rName))
        throw std::invalid_argument(std::string(pOp) + ": no point_data '" + rName + "'");
    const NDArray& a = rMesh.PointData(rName);
    const std::size_t n = rMesh.NumPoints();
    const std::size_t width = n ? a.Size() / n : 0;
    if (a.Size() != n * width || (width != 6 && width != 3) || (Dim == 3 && width != 6))
        throw std::invalid_argument(std::string(pOp) + ": the metric '" + rName +
                                    "' needs 6 components per point (3 in 2-D), the lower "
                                    "triangle m11 m21 m22 m31 m32 m33");
    const detail::DoubleView v(a);
    std::vector<double> out;
    out.reserve(n * (Dim == 3 ? 6 : 3));
    for (std::size_t i = 0; i < n; ++i) {
        double f[6];
        for (std::size_t k = 0; k < width; ++k)
            f[k] = v[i * width + k];
        if (Dim == 3) {
            // file m11 m21 m22 m31 m32 m33 -> API m11 m12 m13 m22 m23 m33
            const double api[6] = {f[0], f[1], f[3], f[2], f[4], f[5]};
            out.insert(out.end(), api, api + 6);
        } else {
            // the xy block of a 3-D metric, or a 2-D one: m11 m12 m22 either way
            out.insert(out.end(), f, f + 3);
        }
    }
    return out;
}

/// Run one MMG entry point under the guard; MMG5_STRONGFAILURE throws.
MmgStatus mmg_run(int Ier, const char* pOp, std::vector<std::string>& rWarnings) {
    if (Ier == MMG5_STRONGFAILURE)
        throw std::runtime_error(std::string(pOp) +
                                 ": MMG failed (MMG5_STRONGFAILURE) and returned no mesh");
    if (Ier == MMG5_LOWFAILURE) {
        rWarnings.push_back(std::string(pOp) +
                            ": MMG returned a valid mesh it could not finish (MMG5_LOWFAILURE)");
        return MmgStatus::LowFailure;
    }
    return MmgStatus::Success;
}

/// What an operation adds to reading MMG's mesh back, by MMG's own references.
struct MmgDecodeExtras {
    /// Elements to leave out (a domain's holes); their orphaned points go too.
    std::vector<MMG5_int> mDropElementRefs;
    /// Facets with this reference become the side region `mInterfaceName`,
    /// on the element whose reference is in `mInsideRefs` when both own it.
    MMG5_int mInterfaceRef = 0;
    std::string mInterfaceName;
    std::vector<MMG5_int> mInsideRefs;
    /// Named cell regions of the elements carrying the given references.
    std::vector<std::pair<std::string, std::vector<MMG5_int> > > mElementGroups;
};

/// Read MMG's mesh back into a meshio++ mesh, with the input's regions.
Mesh mmg_decode(MmgSession& rS, const MmgInput& rIn, const Mesh& rMesh, MmgReport& rReport,
                const char* pOp, const MmgDecodeExtras& rX = {}) {
    MMG5_int np = 0, ntet = 0, nprism = 0, ntri = 0, nquad = 0, nedge = 0;
    MMG5_pMesh m = rS.mMesh;
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            mmg_check(MMG3D_Get_meshSize(m, &np, &ntet, &nprism, &ntri, &nquad, &nedge),
                      "to give the mesh size");
            break;
        case MmgLib::Mmgs:
            mmg_check(MMGS_Get_meshSize(m, &np, &ntri, &nedge), "to give the mesh size");
            break;
        case MmgLib::Mmg2d:
            mmg_check(MMG2D_Get_meshSize(m, &np, &ntri, &nquad, &nedge), "to give the mesh size");
            break;
    }
    const auto unp = static_cast<std::size_t>(np);
    std::vector<double> coords(unp * static_cast<std::size_t>(rIn.mDim));
    std::vector<MMG5_int> vref(unp);
    std::vector<int> corner(unp), req(unp);
    std::vector<MMG5_int> tets(static_cast<std::size_t>(ntet) * 4), tetref(ntet);
    std::vector<MMG5_int> prisms(static_cast<std::size_t>(nprism) * 6), prismref(nprism);
    std::vector<MMG5_int> tris(static_cast<std::size_t>(ntri) * 3), triref(ntri);
    std::vector<MMG5_int> quads(static_cast<std::size_t>(nquad) * 4), quadref(nquad);
    std::vector<MMG5_int> edges(static_cast<std::size_t>(nedge) * 2), edgeref(nedge);
    std::vector<int> elreq(std::max<MMG5_int>({ntet, nprism, ntri, nquad, nedge, 1}));
    std::vector<int> ridge(std::max<MMG5_int>(nedge, 1));
    switch (rS.mLib) {
        case MmgLib::Mmg3d:
            mmg_check(MMG3D_Get_vertices(m, coords.data(), vref.data(), corner.data(), req.data()),
                      "to give the vertices");
            if (ntet)
                mmg_check(MMG3D_Get_tetrahedra(m, tets.data(), tetref.data(), elreq.data()),
                          "to give the tetrahedra");
            if (nprism)
                mmg_check(MMG3D_Get_prisms(m, prisms.data(), prismref.data(), elreq.data()),
                          "to give the prisms");
            if (ntri)
                mmg_check(MMG3D_Get_triangles(m, tris.data(), triref.data(), elreq.data()),
                          "to give the triangles");
            if (nquad)
                mmg_check(MMG3D_Get_quadrilaterals(m, quads.data(), quadref.data(), elreq.data()),
                          "to give the quadrilaterals");
            if (nedge)
                mmg_check(
                    MMG3D_Get_edges(m, edges.data(), edgeref.data(), ridge.data(), elreq.data()),
                    "to give the edges");
            break;
        case MmgLib::Mmgs:
            mmg_check(MMGS_Get_vertices(m, coords.data(), vref.data(), corner.data(), req.data()),
                      "to give the vertices");
            mmg_check(MMGS_Get_triangles(m, tris.data(), triref.data(), elreq.data()),
                      "to give the triangles");
            if (nedge)
                mmg_check(
                    MMGS_Get_edges(m, edges.data(), edgeref.data(), ridge.data(), elreq.data()),
                    "to give the edges");
            break;
        case MmgLib::Mmg2d:
            mmg_check(MMG2D_Get_vertices(m, coords.data(), vref.data(), corner.data(), req.data()),
                      "to give the vertices");
            if (ntri)
                mmg_check(MMG2D_Get_triangles(m, tris.data(), triref.data(), elreq.data()),
                          "to give the triangles");
            if (nquad)
                mmg_check(MMG2D_Get_quadrilaterals(m, quads.data(), quadref.data(), elreq.data()),
                          "to give the quadrilaterals");
            if (nedge)
                mmg_check(
                    MMG2D_Get_edges(m, edges.data(), edgeref.data(), ridge.data(), elreq.data()),
                    "to give the edges");
            break;
    }

    // Leave out the elements of dropped references, then the points no entity uses.
    if (!rX.mDropElementRefs.empty()) {
        const auto dropped = [&](MMG5_int Ref) {
            return std::find(rX.mDropElementRefs.begin(), rX.mDropElementRefs.end(), Ref) !=
                   rX.mDropElementRefs.end();
        };
        const auto filter = [&](std::vector<MMG5_int>& rConn, std::vector<MMG5_int>& rRef,
                                std::size_t K) {
            std::size_t w = 0;
            for (std::size_t i = 0; i < rRef.size(); ++i) {
                if (dropped(rRef[i]))
                    continue;
                for (std::size_t j = 0; j < K; ++j)
                    rConn[w * K + j] = rConn[i * K + j];
                rRef[w++] = rRef[i];
            }
            rRef.resize(w);
            rConn.resize(w * K);
        };
        if (rS.mLib == MmgLib::Mmg3d) {
            filter(tets, tetref, 4);
            filter(prisms, prismref, 6);
        } else {
            filter(tris, triref, 3);
            filter(quads, quadref, 4);
        }
        std::vector<MMG5_int> renumber(unp + 1, 0);
        for (auto* pConn : {&tets, &prisms, &tris, &quads, &edges})
            for (const MMG5_int v : *pConn)
                renumber[static_cast<std::size_t>(v)] = 1;
        MMG5_int next = 0;
        for (std::size_t i = 1; i <= unp; ++i)
            if (renumber[i]) {
                renumber[i] = ++next;
                const std::size_t to = static_cast<std::size_t>(next - 1);
                for (int c = 0; c < rIn.mDim; ++c)
                    coords[to * rIn.mDim + c] = coords[(i - 1) * rIn.mDim + c];
                vref[to] = vref[i - 1];
                corner[to] = corner[i - 1];
                req[to] = req[i - 1];
            }
        for (auto* pConn : {&tets, &prisms, &tris, &quads, &edges})
            for (MMG5_int& v : *pConn)
                v = renumber[static_cast<std::size_t>(v)];
        np = next;
    }
    const auto unp_out = static_cast<std::size_t>(np);

    Mesh out;
    const std::size_t d = rIn.mPointDim;
    NDArray points(DType::Float64, {unp_out, d});
    double* pp = points.As<double>();
    for (std::size_t i = 0; i < unp_out; ++i) {
        for (int c = 0; c < rIn.mDim; ++c)
            pp[i * d + c] = coords[i * rIn.mDim + c];
        if (rIn.mDim == 2 && d == 3)
            pp[i * d + 2] = rIn.mPlaneZ;
    }
    out.AssignPoints(std::move(points));

    // The blocks: elements first, then the facets the input had.
    struct Block {
        const char* mType;
        std::vector<MMG5_int>* pConn;
        std::vector<MMG5_int>* pRef;
        std::size_t mK;
        bool mFacet;
    };
    std::vector<Block> blocks;
    if (rS.mLib == MmgLib::Mmg3d) {
        blocks.push_back({"tetra", &tets, &tetref, 4, false});
        if (nprism)
            blocks.push_back({"wedge", &prisms, &prismref, 6, false});
        blocks.push_back({"triangle", &tris, &triref, 3, true});
        if (nquad)
            blocks.push_back({"quad", &quads, &quadref, 4, true});
        blocks.push_back({"line", &edges, &edgeref, 2, true});
    } else {
        blocks.push_back({"triangle", &tris, &triref, 3, false});
        if (nquad)
            blocks.push_back({"quad", &quads, &quadref, 4, false});
        blocks.push_back({"line", &edges, &edgeref, 2, true});
    }
    // Which facets to keep as cells: the input's own (a reference whose
    // signature is a cell), and MMG's new untagged ones when the input had
    // facets of that kind at all.
    std::vector<std::vector<MMG5_int> > refs_out;
    std::vector<std::vector<MMG5_int> > raw_out;  // MMG's own reference of each cell
    std::vector<bool> facet_out;
    std::size_t num_elements = 0;
    for (const Block& b : blocks) {
        const std::size_t n = b.pRef->size();
        std::vector<std::int64_t> conn;
        std::vector<MMG5_int> refs, raw;
        conn.reserve(n * b.mK);
        for (std::size_t i = 0; i < n; ++i) {
            const MMG5_int ref = (*b.pRef)[i];
            // A reference of the other kind is MMG's own: untagged.
            const MmgSignature* s = rIn.mRefs.Find(ref, b.mFacet);
            if (b.mFacet) {
                const bool lines = std::string(b.mType) == "line";
                const bool had =
                    lines ? !rIn.mEdgeRefs.empty() && rIn.mHadFacetCells : rIn.mHadFacetCells;
                const bool keep = s ? s->mIsCell : had;
                if (!keep)
                    continue;
            }
            for (std::size_t j = 0; j < b.mK; ++j)
                conn.push_back(static_cast<std::int64_t>((*b.pConn)[i * b.mK + j]) - 1);
            refs.push_back(s ? ref : 0);
            raw.push_back(ref);
        }
        if (!b.mFacet)
            num_elements += refs.size();
        if (b.mFacet && refs.empty())
            continue;
        NDArray c(DType::Int64, {refs.size(), b.mK});
        std::copy(conn.begin(), conn.end(), c.As<std::int64_t>());
        out.AddCellBlock(b.mType, std::move(c));
        refs_out.push_back(std::move(refs));
        raw_out.push_back(std::move(raw));
        facet_out.push_back(b.mFacet);
    }

    // medit:ref back, when the input had it.
    if (rIn.mHadPointRef) {
        NDArray pr(DType::Int64, {unp_out});
        for (std::size_t i = 0; i < unp_out; ++i)
            pr.As<std::int64_t>()[i] = vref[i];
        out.AddPointData("medit:ref", std::move(pr));
    }
    // A level set's elements carry MMG's references (its sides, or the
    // materials'), which are what the output's medit:ref then holds.
    if (rIn.mHadCellRef || rIn.mRawElementRefs) {
        std::vector<NDArray> arrays;
        for (std::size_t b = 0; b < refs_out.size(); ++b) {
            const auto& refs = refs_out[b];
            NDArray a(DType::Int64, {refs.size()});
            for (std::size_t i = 0; i < refs.size(); ++i) {
                const MmgSignature* s = rIn.mRefs.Find(refs[i]);
                a.As<std::int64_t>()[i] = rIn.mRawElementRefs && !facet_out[b] ? raw_out[b][i]
                                          : s && s->mHasMeditRef               ? s->mMeditRef
                                                                               : 0;
            }
            arrays.push_back(std::move(a));
        }
        out.AddCellData("medit:ref", std::move(arrays));
    }

    // Cell regions, by signature; side regions through the facets' owners.
    std::vector<std::vector<std::int64_t> > cell_entries(rMesh.NumRegions());
    std::vector<std::vector<std::int64_t> > group_entries(rX.mElementGroups.size());
    std::vector<MMG5_int> element_raw;  // by global cell index, elements only (first)
    std::int64_t global = 0;
    for (std::size_t b = 0; b < refs_out.size(); ++b)
        for (std::size_t i = 0; i < refs_out[b].size(); ++i) {
            if (const MmgSignature* s = rIn.mRefs.Find(refs_out[b][i]))
                for (const std::size_t r : s->mCellRegions)
                    cell_entries[r].push_back(global);
            if (!facet_out[b]) {
                element_raw.push_back(raw_out[b][i]);
                for (std::size_t g = 0; g < rX.mElementGroups.size(); ++g) {
                    const auto& refs = rX.mElementGroups[g].second;
                    if (std::find(refs.begin(), refs.end(), raw_out[b][i]) != refs.end())
                        group_entries[g].push_back(global);
                }
            }
            ++global;
        }
    std::vector<std::int64_t> interface_entries;
    std::vector<std::vector<std::int64_t> > side_entries(rMesh.NumRegions());
    {
        detail::FacetIndexOptions fo;
        fo.mSolidFaces = rS.mLib == MmgLib::Mmg3d;
        fo.mSurfaceEdges = rS.mLib != MmgLib::Mmg3d;
        const detail::FacetIndex index(out, fo);
        const auto facet_sides = [&](const std::vector<MMG5_int>& rConn,
                                     const std::vector<MMG5_int>& rRefs, std::size_t K) {
            for (std::size_t i = 0; i < rRefs.size(); ++i) {
                const MmgSignature* s = rIn.mRefs.Find(rRefs[i], true);
                const bool interface = rX.mInterfaceRef != 0 && rRefs[i] == rX.mInterfaceRef;
                if ((!s || s->mSideRegions.empty()) && !interface)
                    continue;
                std::int64_t corners[4];
                for (std::size_t j = 0; j < K; ++j)
                    corners[j] = static_cast<std::int64_t>(rConn[i * K + j]) - 1;
                const detail::FacetHit* hit = index.Find(corners, K);
                if (!hit)
                    continue;
                if (interface) {
                    // The side of the inside element, when both sides are elements.
                    const auto inside = [&](const detail::FacetOwner& rO) {
                        const auto c = static_cast<std::size_t>(rO.mCell);
                        return c < element_raw.size() &&
                               std::find(rX.mInsideRefs.begin(), rX.mInsideRefs.end(),
                                         element_raw[c]) != rX.mInsideRefs.end();
                    };
                    const detail::FacetOwner& o =
                        hit->mCount > 1 && !inside(hit->mFirst) && inside(hit->mSecond)
                            ? hit->mSecond
                            : hit->mFirst;
                    interface_entries.push_back(o.mCell);
                    interface_entries.push_back(o.mFacet);
                    continue;
                }
                for (const std::size_t r : s->mSideRegions) {
                    side_entries[r].push_back(hit->mFirst.mCell);
                    side_entries[r].push_back(hit->mFirst.mFacet);
                }
            }
        };
        if (rS.mLib == MmgLib::Mmg3d) {
            facet_sides(tris, triref, 3);
            facet_sides(quads, quadref, 4);
        } else {
            facet_sides(edges, edgeref, 2);
        }
    }
    // Point regions: their points survive when required, at their own coordinates.
    std::map<std::vector<double>, std::int64_t> required_out;
    for (std::size_t i = 0; i < unp_out; ++i)
        if (req[i] || corner[i])
            required_out.emplace(std::vector<double>(coords.begin() + i * rIn.mDim,
                                                     coords.begin() + (i + 1) * rIn.mDim),
                                 static_cast<std::int64_t>(i));
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Cell) {
            if (rIn.mRawElementRefs && cell_entries[r].empty() && reg.NumEntries() > 0) {
                rReport.mWarnings.push_back(std::string(pOp) + ": cell region '" + reg.mName +
                                            "' is dropped: the level set gives the elements "
                                            "new references (use materials)");
                continue;
            }
            NDArray e(DType::Int64, {cell_entries[r].size()});
            std::copy(cell_entries[r].begin(), cell_entries[r].end(), e.As<std::int64_t>());
            out.AddRegion(Region(reg.mName, reg.mKind, reg.mDim, reg.mTag, std::move(e)));
            // A domain's boundary lines are also its triangles' sides.
            if (!side_entries[r].empty()) {
                NDArray se(DType::Int64, {side_entries[r].size() / 2, 2});
                std::copy(side_entries[r].begin(), side_entries[r].end(), se.As<std::int64_t>());
                out.AddRegion(
                    Region(reg.mName, RegionKind::Side, reg.mDim, reg.mTag, std::move(se)));
            }
        } else if (reg.mKind == RegionKind::Side) {
            NDArray e(DType::Int64, {side_entries[r].size() / 2, 2});
            std::copy(side_entries[r].begin(), side_entries[r].end(), e.As<std::int64_t>());
            out.AddRegion(Region(reg.mName, reg.mKind, reg.mDim, reg.mTag, std::move(e)));
        } else {
            std::vector<std::int64_t> ids;
            bool all = true;
            const std::int64_t* pe = reg.Entries();
            for (std::size_t k = 0; k < reg.NumEntries(); ++k) {
                const std::size_t p = static_cast<std::size_t>(pe[k]);
                if (p >= rIn.mRequiredPoint.size() || !rIn.mRequiredPoint[p]) {
                    all = false;
                    break;
                }
                std::vector<double> key(rIn.mCoords.begin() + p * rIn.mDim,
                                        rIn.mCoords.begin() + (p + 1) * rIn.mDim);
                auto it = required_out.find(key);
                if (it == required_out.end()) {
                    all = false;
                    break;
                }
                ids.push_back(it->second);
            }
            if (!all) {
                rReport.mWarnings.push_back(std::string(pOp) + ": point region '" + reg.mName +
                                            "' is dropped: its points are not all required "
                                            "(frozen or in a required region)");
                continue;
            }
            NDArray e(DType::Int64, {ids.size()});
            std::copy(ids.begin(), ids.end(), e.As<std::int64_t>());
            out.AddRegion(Region(reg.mName, reg.mKind, reg.mDim, reg.mTag, std::move(e)));
        }
    }
    for (std::size_t g = 0; g < rX.mElementGroups.size(); ++g) {
        NDArray e(DType::Int64, {group_entries[g].size()});
        std::copy(group_entries[g].begin(), group_entries[g].end(), e.As<std::int64_t>());
        out.AddRegion(Region(rX.mElementGroups[g].first, RegionKind::Cell, std::move(e)));
    }
    if (!rX.mInterfaceName.empty()) {
        NDArray e(DType::Int64, {interface_entries.size() / 2, 2});
        std::copy(interface_entries.begin(), interface_entries.end(), e.As<std::int64_t>());
        out.AddRegion(Region(rX.mInterfaceName, RegionKind::Side, std::move(e)));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    rReport.mNumPointsOut = static_cast<std::int64_t>(unp_out);
    rReport.mNumElementsOut = static_cast<std::int64_t>(num_elements);
    const QualityReport q = compute_quality(out);
    for (const auto& [name, summary] : q.mMetrics)
        if (name == "quality:aspect_ratio")
            rReport.mQuality = summary;
    return out;
}

}  // namespace

AdaptResult adapt(const Mesh& rMesh, const AdaptOptions& rOptions) {
    const char* op = "adapt";
    const int given = (rOptions.mHsiz > 0 ? 1 : 0) + (rOptions.mSize.empty() ? 0 : 1) +
                      (rOptions.mMetric.empty() ? 0 : 1);
    if (given > 1)
        throw std::invalid_argument("adapt: give one of hsiz, size and metric");
    AdaptResult result;
    MmgReport& report = result.mReport;
    MmgInput in = mmg_encode(rMesh, rOptions.mParameters, op, {rOptions.mSize, rOptions.mMetric},
                             report.mWarnings);
    report.mLibrary = mmg_lib_name(in.mLib);
    report.mNumPointsIn = static_cast<std::int64_t>(rMesh.NumPoints());
    report.mNumElementsIn = static_cast<std::int64_t>(in.mNumElements);
    std::vector<double> sol;
    int sol_type = 0;
    if (!rOptions.mSize.empty()) {
        sol = mmg_scalar_field(rMesh, rOptions.mSize, op);
        if (std::any_of(sol.begin(), sol.end(), [](double V) { return !(V > 0.0); }))
            throw std::invalid_argument("adapt: the size field must be positive everywhere");
        sol_type = MMG5_Scalar;
    } else if (!rOptions.mMetric.empty()) {
        sol = mmg_metric_field(rMesh, rOptions.mMetric, in.mDim, op);
        sol_type = MMG5_Tensor;
    }

    {
        detail::MmgCallGuard guard;
        MmgSession s(in.mLib, false, false);
        mmg_load(s, in);
        mmg_parameters(s, in, rMesh, rOptions.mParameters, op);
        if (rOptions.mHsiz > 0)
            mmg_check(mmg_set_d(s, s.mMet, MMG3D_DPARAM_hsiz, MMGS_DPARAM_hsiz, MMG2D_DPARAM_hsiz,
                                rOptions.mHsiz),
                      "the parameter hsiz");
        if (sol_type)
            mmg_vertex_sol(s, s.mMet, sol_type, sol);
        int ier = MMG5_STRONGFAILURE;
        switch (in.mLib) {
            case MmgLib::Mmg3d:
                ier = MMG3D_mmg3dlib(s.mMesh, s.mMet);
                break;
            case MmgLib::Mmgs:
                ier = MMGS_mmgslib(s.mMesh, s.mMet);
                break;
            case MmgLib::Mmg2d:
                ier = MMG2D_mmg2dlib(s.mMesh, s.mMet);
                break;
        }
        report.mStatus = mmg_run(ier, op, report.mWarnings);
        result.mMesh = mmg_decode(s, in, rMesh, report, op);
    }
    for (const std::string& w : report.mWarnings)
        log::warn("{}", w);
    return result;
}

DiscretizeLevelsetResult discretize_levelset(const Mesh& rMesh,
                                             const DiscretizeLevelsetOptions& rOptions) {
    const char* op = "discretize_levelset";
    if (rOptions.mLevelSet.empty())
        throw std::invalid_argument("discretize_levelset: name the level set, a point_data array");
    DiscretizeLevelsetResult result;
    MmgReport& report = result.mReport;
    MmgInput in = mmg_encode(rMesh, rOptions.mParameters, op, {rOptions.mLevelSet},
                             report.mWarnings, MmgEncodeMode::LevelSet);
    report.mLibrary = mmg_lib_name(in.mLib);
    report.mNumPointsIn = static_cast<std::int64_t>(rMesh.NumPoints());
    report.mNumElementsIn = static_cast<std::int64_t>(in.mNumElements);
    std::vector<double> ls = mmg_scalar_field(rMesh, rOptions.mLevelSet, op);

    MmgDecodeExtras x;
    x.mInterfaceRef = MG_ISO;
    x.mInterfaceName = "interface";
    std::vector<MMG5_int> inside = {MG_MINUS}, outside = {MG_PLUS};
    for (const MmgMaterial& m : rOptions.mMaterials) {
        inside.push_back(static_cast<MMG5_int>(m.mRefInside));
        outside.push_back(static_cast<MMG5_int>(m.mRefOutside));
    }
    x.mInsideRefs = inside;
    x.mElementGroups = {{"inside", inside}, {"outside", outside}};
    {
        detail::MmgCallGuard guard;
        MmgSession s(in.mLib, true, false);
        mmg_load(s, in);
        mmg_parameters(s, in, rMesh, rOptions.mParameters, op);
        mmg_check(mmg_set_i(s, s.mMet, MMG3D_IPARAM_iso, MMGS_IPARAM_iso, MMG2D_IPARAM_iso, 1),
                  "the parameter iso");
        mmg_check(mmg_set_d(s, s.mMet, MMG3D_DPARAM_ls, MMGS_DPARAM_ls, MMG2D_DPARAM_ls,
                            rOptions.mIsovalue),
                  "the isovalue");
        if (rOptions.mRemoveSmallComponents > 0)
            mmg_check(mmg_set_d(s, s.mMet, MMG3D_DPARAM_rmc, MMGS_DPARAM_rmc, MMG2D_DPARAM_rmc,
                                rOptions.mRemoveSmallComponents),
                      "the parameter rmc");
        if (!rOptions.mMaterials.empty()) {
            mmg_check(mmg_set_i(s, s.mMet, MMG3D_IPARAM_numberOfMat, MMGS_IPARAM_numberOfMat,
                                MMG2D_IPARAM_numberOfMat,
                                static_cast<MMG5_int>(rOptions.mMaterials.size())),
                      "the number of materials");
            for (const MmgMaterial& m : rOptions.mMaterials) {
                const int split = m.mSplit ? MMG5_MMAT_Split : MMG5_MMAT_NoSplit;
                const auto ref = static_cast<MMG5_int>(m.mRef);
                const auto rin = static_cast<MMG5_int>(m.mRefInside);
                const auto rout = static_cast<MMG5_int>(m.mRefOutside);
                int ok = 0;
                switch (in.mLib) {
                    case MmgLib::Mmg3d:
                        ok = MMG3D_Set_multiMat(s.mMesh, s.mLs, ref, split, rin, rout);
                        break;
                    case MmgLib::Mmgs:
                        ok = MMGS_Set_multiMat(s.mMesh, s.mLs, ref, split, rin, rout);
                        break;
                    case MmgLib::Mmg2d:
                        ok = MMG2D_Set_multiMat(s.mMesh, s.mLs, ref, split, rin, rout);
                        break;
                }
                mmg_check(ok, "a material");
            }
        }
        mmg_vertex_sol(s, s.mLs, MMG5_Scalar, ls);
        int ier = MMG5_STRONGFAILURE;
        switch (in.mLib) {
            case MmgLib::Mmg3d:
                ier = MMG3D_mmg3dls(s.mMesh, s.mLs, s.mMet);
                break;
            case MmgLib::Mmgs:
                ier = MMGS_mmgsls(s.mMesh, s.mLs, s.mMet);
                break;
            case MmgLib::Mmg2d:
                ier = MMG2D_mmg2dls(s.mMesh, s.mLs, s.mMet);
                break;
        }
        report.mStatus = mmg_run(ier, op, report.mWarnings);
        result.mMesh = mmg_decode(s, in, rMesh, report, op, x);
    }
    for (const std::string& w : report.mWarnings)
        log::warn("{}", w);
    return result;
}

MeshDomainResult mesh_domain(const Mesh& rBoundary, const MeshDomainOptions& rOptions) {
    const char* op = "mesh_domain";
    if (rOptions.mHsiz > 0 && !rOptions.mSize.empty())
        throw std::invalid_argument("mesh_domain: give one of hsiz and size");
    if (rOptions.mHoles.size() % 2 != 0)
        throw std::invalid_argument("mesh_domain: holes are x, y pairs");
    MeshDomainResult result;
    MmgReport& report = result.mReport;
    MmgInput in = mmg_encode(rBoundary, rOptions.mParameters, op, {rOptions.mSize},
                             report.mWarnings, MmgEncodeMode::Domain);
    report.mLibrary = mmg_lib_name(in.mLib);
    report.mNumPointsIn = static_cast<std::int64_t>(rBoundary.NumPoints());
    std::vector<double> size;
    if (!rOptions.mSize.empty()) {
        size = mmg_scalar_field(rBoundary, rOptions.mSize, op);
        if (std::any_of(size.begin(), size.end(), [](double V) { return !(V > 0.0); }))
            throw std::invalid_argument("mesh_domain: the size field must be positive everywhere");
    }
    {
        detail::MmgCallGuard guard;
        MmgSession s(in.mLib, false, false);
        mmg_load(s, in);
        mmg_parameters(s, in, rBoundary, rOptions.mParameters, op);
        if (rOptions.mHsiz > 0)
            mmg_check(MMG2D_Set_dparameter(s.mMesh, s.mMet, MMG2D_DPARAM_hsiz, rOptions.mHsiz),
                      "the parameter hsiz");
        if (!size.empty())
            mmg_vertex_sol(s, s.mMet, MMG5_Scalar, size);
        report.mStatus = mmg_run(MMG2D_mmg2dmesh(s.mMesh, s.mMet), op, report.mWarnings);

        // MMG gives each enclosed region its own reference and has no holes:
        // find the regions the hole points are in.
        MMG5_int np = 0, nt = 0, nquad = 0, na = 0;
        mmg_check(MMG2D_Get_meshSize(s.mMesh, &np, &nt, &nquad, &na), "to give the mesh size");
        std::vector<double> xy(static_cast<std::size_t>(np) * 2);
        std::vector<MMG5_int> vref(static_cast<std::size_t>(np));
        std::vector<int> flags(static_cast<std::size_t>(std::max<MMG5_int>(np, nt) + 1));
        std::vector<int> flags2(flags.size());
        mmg_check(MMG2D_Get_vertices(s.mMesh, xy.data(), vref.data(), flags.data(), flags2.data()),
                  "to give the vertices");
        std::vector<MMG5_int> tri(static_cast<std::size_t>(nt) * 3),
            tref(static_cast<std::size_t>(nt));
        if (nt)
            mmg_check(MMG2D_Get_triangles(s.mMesh, tri.data(), tref.data(), flags.data()),
                      "to give the triangles");
        MmgDecodeExtras x;
        for (std::size_t h = 0; h + 1 < rOptions.mHoles.size(); h += 2) {
            const double px = rOptions.mHoles[h], py = rOptions.mHoles[h + 1];
            MMG5_int found = 0;
            for (std::size_t t = 0; t < static_cast<std::size_t>(nt) && !found; ++t) {
                const double* a = &xy[static_cast<std::size_t>(tri[t * 3] - 1) * 2];
                const double* b = &xy[static_cast<std::size_t>(tri[t * 3 + 1] - 1) * 2];
                const double* c = &xy[static_cast<std::size_t>(tri[t * 3 + 2] - 1) * 2];
                const auto side = [](const double* P, const double* Q, double X, double Y) {
                    return (Q[0] - P[0]) * (Y - P[1]) - (Q[1] - P[1]) * (X - P[0]);
                };
                const double d0 = side(a, b, px, py), d1 = side(b, c, px, py),
                             d2 = side(c, a, px, py);
                if ((d0 >= 0 && d1 >= 0 && d2 >= 0) || (d0 <= 0 && d1 <= 0 && d2 <= 0))
                    found = tref[t];
            }
            if (!found)
                throw std::invalid_argument("mesh_domain: the hole point (" + std::to_string(px) +
                                            ", " + std::to_string(py) + ") is outside the domain");
            if (std::find(x.mDropElementRefs.begin(), x.mDropElementRefs.end(), found) ==
                x.mDropElementRefs.end())
                x.mDropElementRefs.push_back(found);
        }
        std::vector<MMG5_int> kept;
        for (const MMG5_int ref : tref)
            if (std::find(kept.begin(), kept.end(), ref) == kept.end() &&
                std::find(x.mDropElementRefs.begin(), x.mDropElementRefs.end(), ref) ==
                    x.mDropElementRefs.end())
                kept.push_back(ref);
        std::sort(kept.begin(), kept.end());
        for (const MMG5_int ref : kept)
            x.mElementGroups.push_back({"subdomain:" + std::to_string(ref), {ref}});
        result.mNumSubdomains = static_cast<std::int64_t>(kept.size());
        result.mNumHoles = static_cast<std::int64_t>(x.mDropElementRefs.size());
        result.mMesh = mmg_decode(s, in, rBoundary, report, op, x);
    }
    for (const std::string& w : report.mWarnings)
        log::warn("{}", w);
    return result;
}

AdaptResult move(const Mesh& rMesh, const MoveOptions& rOptions) {
    const char* op = "move";
#ifndef MESHIOPLUSPLUS_HAS_MMG_ELAS
    (void)rMesh;
    (void)rOptions;
    throw Unsupported(
        "move: the linked MMG was built without ELAS, which its Lagrangian mode needs: "
        "rebuild MMG with build/build-mmg.sh --with-elas (GPL-3.0; see doc/mmg.md)");
#else
    if (rOptions.mDisplacement.empty())
        throw std::invalid_argument("move: name the displacement, a point_data vector");
    if (rOptions.mMoving.empty())
        throw std::invalid_argument("move: name the regions whose boundary moves");
    if (rOptions.mLagrangian < 0 || rOptions.mLagrangian > 2)
        throw std::invalid_argument("move: lagrangian must be 0, 1 or 2");
    AdaptResult result;
    MmgReport& report = result.mReport;
    MmgInput in = mmg_encode(rMesh, rOptions.mParameters, op, {rOptions.mDisplacement},
                             report.mWarnings, MmgEncodeMode::Move);
    if (in.mLib == MmgLib::Mmgs)
        throw std::invalid_argument(
            "move: MMG moves tetrahedra (mmg3d) and planar triangles (mmg2d), not a surface");
    std::vector<std::size_t> moving;
    for (const std::string& rName : rOptions.mMoving) {
        bool found = false;
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
            if (rMesh.Region(r).mName == rName && rMesh.Region(r).mKind != RegionKind::Point) {
                moving.push_back(r);
                found = true;
            }
        if (!found)
            throw std::invalid_argument("move: no cell or side region named '" + rName + "'");
    }
    report.mLibrary = mmg_lib_name(in.mLib);
    report.mNumPointsIn = static_cast<std::int64_t>(rMesh.NumPoints());
    report.mNumElementsIn = static_cast<std::int64_t>(in.mNumElements);

    // The moving facets carry MMG5_DISPREF (10, private to MMG:
    // mmgcommon_private.h); their own reference (merged when they had
    // several) is put back on MMG's output.
    constexpr MMG5_int kMmgDispRef = 10;
    std::vector<MMG5_int>& refs = in.mLib == MmgLib::Mmg3d ? in.mTriRefs : in.mEdgeRefs;
    MmgSignature merged;
    merged.mFacet = true;
    std::vector<MMG5_int> originals;
    std::size_t num_moving = 0;
    for (MMG5_int& rRef : refs) {
        const MmgSignature* sig = in.mRefs.Find(rRef, true);
        if (!sig)
            continue;
        const bool hit = std::any_of(moving.begin(), moving.end(), [&](std::size_t R) {
            return std::find(sig->mCellRegions.begin(), sig->mCellRegions.end(), R) !=
                       sig->mCellRegions.end() ||
                   std::find(sig->mSideRegions.begin(), sig->mSideRegions.end(), R) !=
                       sig->mSideRegions.end();
        });
        if (!hit)
            continue;
        if (std::find(originals.begin(), originals.end(), rRef) == originals.end()) {
            originals.push_back(rRef);
            merged.mIsCell = merged.mIsCell || sig->mIsCell;
            if (originals.size() == 1) {
                merged.mHasMeditRef = sig->mHasMeditRef;
                merged.mMeditRef = sig->mMeditRef;
            } else if (!sig->mHasMeditRef || sig->mMeditRef != merged.mMeditRef) {
                merged.mHasMeditRef = false;
            }
            for (const std::size_t r : sig->mCellRegions)
                mmg_add_unique(merged.mCellRegions, r);
            for (const std::size_t r : sig->mSideRegions)
                mmg_add_unique(merged.mSideRegions, r);
        }
        rRef = kMmgDispRef;
        ++num_moving;
    }
    if (num_moving == 0)
        throw std::invalid_argument(
            "move: the moving regions hold no boundary facet (a side region, or facet cells)");
    const MMG5_int restored = originals.size() == 1 ? originals[0] : in.mRefs.Ref(merged);
    if (originals.size() > 1)
        report.mWarnings.push_back(
            "move: the moving facets had several references; they come back with the union of "
            "their regions");

    // The displacement, dim components per point (xy of a 3-component one for mmg2d).
    if (!rMesh.HasPointData(rOptions.mDisplacement))
        throw std::invalid_argument("move: no point_data '" + rOptions.mDisplacement + "'");
    const NDArray& d = rMesh.PointData(rOptions.mDisplacement);
    const std::size_t np = rMesh.NumPoints();
    const std::size_t width = np ? d.Size() / np : 0;
    const auto dim = static_cast<std::size_t>(in.mDim);
    if (d.Size() != np * width || width < dim || width > 3)
        throw std::invalid_argument("move: the displacement must be a point_data vector of " +
                                    std::to_string(dim) + " components (shape (n, " +
                                    std::to_string(dim) + "))");
    const detail::DoubleView flat(d);
    std::vector<double> disp(np * dim);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < dim; ++k)
            disp[i * dim + k] = flat[i * width + k];

    {
        detail::MmgCallGuard guard;
        // ELAS prints whatever the verbosity (detail/mmg_guard.hpp).
        detail::MmgStdoutSilencer quiet(rOptions.mParameters.mVerbose < 0);
        MmgSession s(in.mLib, false, true);
        mmg_load(s, in);
        mmg_parameters(s, in, rMesh, rOptions.mParameters, op);
        mmg_check(mmg_set_i(s, s.mMet, MMG3D_IPARAM_lag, -1, MMG2D_IPARAM_lag,
                            static_cast<MMG5_int>(rOptions.mLagrangian)),
                  "the parameter lag");
        mmg_vertex_sol(s, s.mDisp, MMG5_Vector, disp);
        const int ier = in.mLib == MmgLib::Mmg3d ? MMG3D_mmg3dmov(s.mMesh, s.mMet, s.mDisp)
                                                 : MMG2D_mmg2dmov(s.mMesh, s.mMet, s.mDisp);
        report.mStatus = mmg_run(ier, op, report.mWarnings);
        MMG5_pMesh m = s.mMesh;
        if (in.mLib == MmgLib::Mmg3d) {
            for (MMG5_int k = 1; k <= m->nt; ++k)
                if (m->tria[k].ref == kMmgDispRef)
                    m->tria[k].ref = restored;
        } else {
            for (MMG5_int k = 1; k <= m->na; ++k)
                if (m->edge[k].ref == kMmgDispRef)
                    m->edge[k].ref = restored;
        }
        result.mMesh = mmg_decode(s, in, rMesh, report, op);
    }
    for (const std::string& w : report.mWarnings)
        log::warn("{}", w);
    return result;
#endif
}

#endif  // MESHIOPLUSPLUS_HAS_MMG

}  // namespace meshioplusplus
