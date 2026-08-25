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
// See gid.hpp for the format's full documentation.

// System includes
#include <stdexcept>

// Project includes
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/exceptions.hpp"

#include "gid_common.hpp"

namespace meshioplusplus {

// Unconditional (needed whether or not gidpost itself is compiled in), the
// smooth_method_from_name/partition_method_from_name precedent: the
// flat-binding/CLI/Python spelling of GidMode.
GidMode gid_mode_from_name(const std::string& rName) {
    if (rName == "auto")
        return GidMode::Auto;
    if (rName == "ascii")
        return GidMode::Ascii;
    if (rName == "binary")
        return GidMode::Binary;
    if (rName == "hdf5")
        return GidMode::Hdf5;
    throw std::invalid_argument("GiD: unknown mode '" + rName +
                                "' (expected auto, ascii, binary, or hdf5)");
}

}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_HAS_GIDPOST

// System includes
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"

extern "C" {
#include "gidpost.h"
}

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Cell-type mapping.
//
// GiD has exactly ten element types (GiD_ElementType); higher-order variants
// share a type with a larger node count. Every entry below is IDENTITY (no
// node permutation) -- independently cross-checked against the real,
// production GiD writer in Kratos Multiphysics (kratos/includes/
// gid_mesh_container.h), whose only reordering is for Hexahedra20, and that
// reorder turns out to exist purely because KRATOS's own internal hexahedra20
// node order differs from GiD's (Kratos: corners, bottom-ring, verticals,
// top-ring; GiD: corners, bottom-ring, top-ring, verticals) -- and GiD's own
// convention is, edge for edge, IDENTICAL to meshio++'s own hexahedron20
// table (detail/cell_subdivision.cpp's cell_refine_edges(Hexahedron)):
// bottom ring {0,1}{1,2}{2,3}{3,0}, top ring {4,5}{5,6}{6,7}{7,4}, verticals
// {0,4}{1,5}{2,6}{3,7}. Kratos's own geometry classes (triangle_2d_6.h,
// tetrahedra_3d_10.h, quadrilateral_2d_8.h) write NO reorder at all for
// triangle6/tetra10/quad8, and their own edge tables match meshio++'s
// (detail/cell_edges.cpp) edge-for-edge too -- so those three, and quad9 by
// the same near-universal "trailing node = face centre" convention, are
// identity as well. Pinned independently of this derivation by the
// tests/cpp/test_gid.cpp GidOrdering suite, which checks the RAW WRITTEN
// FILE against GiD's own geometry, never a round trip (there is no reader).
//
// UNRESOLVED, and deliberately recorded rather than quietly settled:
// CIMNE's own published figure for the 20-node hexahedron (the `hexa20.gif`
// in the GiD reference manual's postprocess-format page) numbers the
// mid-edge nodes bottom-ring, VERTICALS, top-ring -- i.e. exactly Kratos's
// INTERNAL order, the order Kratos then permutes away from. Taken at face
// value it says the identity mapping used here is wrong. It is not followed,
// for two reasons: the figure is from the GiD 6-era manual and CIMNE's
// current grammar dropped the mid-edge figures entirely, saying only
// "hierarchical order ... vertex nodes first, then the middle ones"; and
// Kratos's permutation is a production code path exercised against real GiD
// for years and labelled a "workaround", i.e. added in response to an
// observed problem. Documentary evidence loses to that. Settling it needs an
// external oracle nobody here has -- a hexahedron20 file written by GiD
// itself, or GiD rendering ours -- so the risk is stated in
// doc/formats/gid.md instead of being hidden behind a confident comment.
// NOTE the GidOrdering tests cannot adjudicate this: they pin that no
// permutation is applied, not that none is needed.
//
// Anything not in this table -- hexahedron27/wedge15/pyramid13 (orderings
// not independently verified above), polygon/polyhedron (GiD has no such
// type), every VTK-Lagrange/higher-degree type -- throws by name rather than
// guessing.
struct GidTypeEntry {
    GiD_ElementType mType;
    int mNumNodes;
};

const std::unordered_map<std::string, GidTypeEntry>& gid_type_table() {
    static const std::unordered_map<std::string, GidTypeEntry> table = {
        {"vertex", {GiD_Point, 1}},
        {"line", {GiD_Linear, 2}},
        {"line3", {GiD_Linear, 3}},
        {"triangle", {GiD_Triangle, 3}},
        {"triangle6", {GiD_Triangle, 6}},
        {"quad", {GiD_Quadrilateral, 4}},
        {"quad8", {GiD_Quadrilateral, 8}},
        {"quad9", {GiD_Quadrilateral, 9}},
        {"tetra", {GiD_Tetrahedra, 4}},
        {"tetra10", {GiD_Tetrahedra, 10}},
        {"hexahedron", {GiD_Hexahedra, 8}},
        {"hexahedron20", {GiD_Hexahedra, 20}},
        {"wedge", {GiD_Prism, 6}},
        {"pyramid", {GiD_Pyramid, 5}},
    };
    return table;
}

bool gid_is_int_dtype(DType t) {
    return t == DType::Int8 || t == DType::Int16 || t == DType::Int32 || t == DType::Int64 ||
           t == DType::UInt8 || t == DType::UInt16 || t == DType::UInt32 || t == DType::UInt64;
}

// gidpost validates its own state machine with assert() throughout, and
// GiD_PostInit() prints a "Debug version" banner to stdout in a non-NDEBUG
// build -- both compiled out here via NDEBUG (CMakeLists.txt's
// _mio_gidpost_src_defs), which is what makes the return-code checks below
// the ONLY error-reporting path: an abort() inside the Python extension or
// the WASM module would kill the host instead of raising, contradicting this
// project's errors-are-exceptions contract.
void gid_ensure_init() {
    static const bool once = [] {
        GiD_PostInit();
        return true;
    }();
    (void)once;
    // GiD_PostDone() is deliberately never called: it tears down gidpost's
    // process-global file-handle table, which a later write_gid() call in
    // the same process would need.
}

void gid_check(int rc, const std::string& rWhat) {
    if (rc != 0)
        throw WriteError("GiD: " + rWhat + " failed (gidpost status " + std::to_string(rc) + ")");
}

// gid_has_suffix / gid_resolve_mode / gid_ascii_paths now live in the
// format-private formats/gid_common.hpp, so the reader -- which must NOT be
// behind this file's MESHIOPLUSPLUS_HAS_GIDPOST guard -- can share them.
using gid_detail::gid_ascii_paths;
using gid_detail::gid_resolve_mode;

GiD_PostMode gid_post_mode(GidMode mode) {
    switch (mode) {
        case GidMode::Binary:
            return GiD_PostBinary;
        case GidMode::Hdf5:
            return GiD_PostHDF5;
        default:
            return GiD_PostAscii;
    }
}

// Provenance: SlotTier::Block, rendered as one GiD "user attribute" per line
// (gidpost's own source documents these as rendering `# Name: value` in
// ASCII/raw-binary and an HDF5 attribute otherwise) -- "meshio++" for line 0
// (kProvenanceTag), "meshio++_<n>" for any further line under an open
// BestEffort/Required scope.
void gid_write_provenance_mesh(GiD_FILE fd) {
    const std::vector<std::string> lines = detail::provenance_lines(detail::SlotTier::Block);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string name = i == 0 ? "meshio++" : "meshio++_" + std::to_string(i);
        GiD_fWriteMeshUserAttribute(fd, name.c_str(), lines[i].c_str());
    }
}

void gid_write_provenance_result(GiD_FILE fd) {
    const std::vector<std::string> lines = detail::provenance_lines(detail::SlotTier::Block);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string name = i == 0 ? "meshio++" : "meshio++_" + std::to_string(i);
        GiD_fWriteResultUserAttribute(fd, name.c_str(), lines[i].c_str());
    }
}

// ---------------------------------------------------------------------------
// Geometry.
// ---------------------------------------------------------------------------

void gid_write_geometry(GiD_FILE fd, const Mesh& rMesh,
                        const std::vector<const GidTypeEntry*>& rEntries,
                        const std::vector<std::string>& rMeshNames, const std::string* pMatName) {
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    std::vector<double> xyz(np * 3);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t c = 0; c < 3; ++c)
            xyz[i * 3 + c] = c < dim ? detail::read_double(points, i * dim + c) : 0.0;

    std::int64_t elem_base = 1;  // next 1-based global element id, across every block
    bool wrote_coords = false;

    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const GidTypeEntry* entry = rEntries[bi];
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t ne = cb.NumCells();
        const NDArray& conn = cb.Conn();

        gid_check(GiD_fBeginMesh(fd, rMeshNames[bi].c_str(), GiD_3D, entry->mType,
                                 static_cast<int>(npc)),
                  "BeginMesh('" + rMeshNames[bi] + "')");
        if (bi == 0)
            gid_write_provenance_mesh(fd);

        // Only the first mesh carries the shared node table; gidpost's own
        // state machine requires every subsequent mesh to still open/close an
        // (empty) coordinates block, or it rejects the file.
        if (!wrote_coords) {
            gid_check(GiD_fWriteCoordinatesBlock(fd, static_cast<int>(np), xyz.data()),
                      "WriteCoordinatesBlock");
            wrote_coords = true;
        } else {
            gid_check(GiD_fBeginCoordinates(fd), "BeginCoordinates");
            gid_check(GiD_fEndCoordinates(fd), "EndCoordinates");
        }

        // GiD_fWriteElementsId(Mat)Block already wraps BeginElements()/
        // EndElements() internally (gidpost's own header documents this) --
        // an explicit pair here would double-nest the section.
        std::vector<int> ids(ne);
        std::vector<int> flat_conn(ne * npc);
        for (std::size_t r = 0; r < ne; ++r) {
            ids[r] = static_cast<int>(elem_base + static_cast<std::int64_t>(r));
            for (std::size_t j = 0; j < npc; ++j)
                flat_conn[r * npc + j] = static_cast<int>(detail::read_int(conn, r * npc + j)) + 1;
        }

        if (pMatName != nullptr) {
            const NDArray& mat_arr = rMesh.CellData(*pMatName, bi);
            std::vector<int> mat(ne);
            for (std::size_t r = 0; r < ne; ++r)
                mat[r] = static_cast<int>(detail::read_int(mat_arr, r));
            gid_check(GiD_fWriteElementsIdMatBlock(fd, static_cast<int>(ne), ids.data(),
                                                   flat_conn.data(), mat.data()),
                      "WriteElementsIdMatBlock('" + rMeshNames[bi] + "')");
        } else {
            gid_check(
                GiD_fWriteElementsIdBlock(fd, static_cast<int>(ne), ids.data(), flat_conn.data()),
                "WriteElementsIdBlock('" + rMeshNames[bi] + "')");
        }

        gid_check(GiD_fEndMesh(fd), "EndMesh");

        elem_base += static_cast<std::int64_t>(ne);
    }
}

// ---------------------------------------------------------------------------
// Results.
// ---------------------------------------------------------------------------

// GiD_ResultType's valid component counts are irregular (Scalar 1; Vector
// 2/3/4; Matrix 3/6; MainMatrix 12; ...) and gidpost does not validate an
// unsupported count itself -- it silently emits a malformed file. Only 1
// (Scalar) and 2/3 (Vector) are mapped directly; anything else splits into
// that many named scalars (see gid.hpp's k-component rule for why a
// 6-component array is deliberately not mapped to GiD_Matrix).
void gid_write_result_array(GiD_FILE fd, const NDArray& rArr, std::size_t rows,
                            const std::vector<int>& rIds, const std::string& rName,
                            const std::string& rAnalysis, double step, GiD_ResultLocation loc,
                            const std::string& rGaussName) {
    const std::size_t k = rArr.Shape().size() >= 2 ? rArr.Shape()[1] : 1;
    const char* gauss = rGaussName.empty() ? nullptr : rGaussName.c_str();

    if (k == 1 || k == 2 || k == 3) {
        const GiD_ResultType rtype = k == 1 ? GiD_Scalar : GiD_Vector;
        std::vector<double> vals(rows * k);
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c < k; ++c)
                vals[r * k + c] = detail::read_double(rArr, r * k + c);
        gid_check(GiD_fWriteResultBlock(fd, rName.c_str(), rAnalysis.c_str(), step, rtype, loc,
                                        gauss, nullptr, 0, nullptr, nullptr,
                                        static_cast<int>(rows), rIds.data(), static_cast<int>(k),
                                        vals.data()),
                  "WriteResultBlock('" + rName + "')");
        return;
    }

    detail::provenance_note("result-split", rName + ": k=" + std::to_string(k) +
                                                " has no GiD result type; written as " +
                                                std::to_string(k) + " scalars");
    std::vector<double> col(rows);
    for (std::size_t c = 0; c < k; ++c) {
        for (std::size_t r = 0; r < rows; ++r)
            col[r] = detail::read_double(rArr, r * k + c);
        const std::string cname = rName + "_" + std::to_string(c + 1);
        gid_check(GiD_fWriteResultBlock(fd, cname.c_str(), rAnalysis.c_str(), step, GiD_Scalar, loc,
                                        gauss, nullptr, 0, nullptr, nullptr,
                                        static_cast<int>(rows), rIds.data(), 1, col.data()),
                  "WriteResultBlock('" + cname + "')");
    }
}

void gid_write_point_data(GiD_FILE fd, const Mesh& rMesh, const std::string& rAnalysis,
                          double step) {
    const std::size_t np = rMesh.NumPoints();
    std::vector<int> ids(np);
    for (std::size_t i = 0; i < np; ++i)
        ids[i] = static_cast<int>(i) + 1;

    for (const auto& name : rMesh.PointDataNames())
        gid_write_result_array(fd, rMesh.PointData(name), np, ids, name, rAnalysis, step,
                               GiD_OnNodes, "");
}

// GiD_ResultLocation has no "on cells" concept, so cell_data is written
// GiD_OnGaussPoints against a synthetic one-Gauss-point set declared once per
// block ("gp_<mesh_name>") -- the standard GiD idiom for a per-element field,
// and how Kratos's own GiD writer represents one.
void gid_write_cell_data(GiD_FILE fd, const Mesh& rMesh,
                         const std::vector<const GidTypeEntry*>& rEntries,
                         const std::vector<std::string>& rMeshNames,
                         const std::vector<std::int64_t>& rElemBase, const std::string* pMatName,
                         const std::string& rAnalysis, double step) {
    const std::size_t nb = rMesh.NumCellBlocks();

    bool any_cell_data = false;
    for (const auto& name : rMesh.CellDataNames())
        if (pMatName == nullptr || name != *pMatName) {
            any_cell_data = true;
            break;
        }
    if (!any_cell_data)
        return;

    // One Gauss-point set per block, declared once regardless of how many
    // arrays reference it.
    for (std::size_t bi = 0; bi < nb; ++bi) {
        const auto cb = rMesh.Cells(bi);
        const std::string gauss_name = "gp_" + rMeshNames[bi];
        gid_check(GiD_fBeginGaussPoint(fd, gauss_name.c_str(), rEntries[bi]->mType,
                                       rMeshNames[bi].c_str(), /*GP_number=*/1,
                                       /*NodesIncluded=*/0, /*InternalCoord=*/1),
                  "BeginGaussPoint('" + gauss_name + "')");
        gid_check(GiD_fEndGaussPoint(fd), "EndGaussPoint('" + gauss_name + "')");
        (void)cb;
    }

    for (const auto& name : rMesh.CellDataNames()) {
        if (pMatName != nullptr && name == *pMatName)
            continue;  // already written as the geometry file's material column
        if (rMesh.CellDataNumBlocks(name) != nb) {
            log::warn("gid: cell_data '{}' has {} blocks, mesh has {} -- skipped", name,
                     rMesh.CellDataNumBlocks(name), nb);
            continue;
        }
        for (std::size_t bi = 0; bi < nb; ++bi) {
            const auto cb = rMesh.Cells(bi);
            const std::size_t ne = cb.NumCells();
            if (ne == 0)
                continue;
            std::vector<int> ids(ne);
            for (std::size_t r = 0; r < ne; ++r)
                ids[r] = static_cast<int>(rElemBase[bi] + static_cast<std::int64_t>(r));
            const std::string gauss_name = "gp_" + rMeshNames[bi];
            gid_write_result_array(fd, rMesh.CellData(name, bi), ne, ids, name, rAnalysis, step,
                                   GiD_OnGaussPoints, gauss_name);
        }
    }
}

}  // namespace

bool gid_available(GidMode mode) {
    if (mode == GidMode::Hdf5) {
#ifdef MESHIOPLUSPLUS_HAS_GIDPOST_HDF5
        return true;
#else
        return false;
#endif
    }
    return true;
}

std::string gid_build_option(GidMode mode) {
    if (mode == GidMode::Hdf5)
        return "-DMESHIOPLUSPLUS_WITH_GIDPOST=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON "
               "-DMESHIOPLUSPLUS_WITH_HDF5=ON";
    return "-DMESHIOPLUSPLUS_WITH_GIDPOST=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON";
}

void write_gid(const std::string& rPath, const Mesh& rMesh, GidMode mode,
              const std::string& rAnalysisName, double stepValue) {
    const GidMode resolved = gid_resolve_mode(rPath, mode);
    if (!gid_available(resolved))
        throw WriteError("meshio++: the 'gid' HDF5 flavour needs a build with " +
                         gid_build_option(resolved));

    if (rMesh.PointDim() > 3)
        throw WriteError("GiD: points must have at most three components");
    constexpr std::size_t kIntMax = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (rMesh.NumPoints() > kIntMax)
        throw WriteError("GiD: mesh has too many points for gidpost's 32-bit API");

    const std::size_t nb = rMesh.NumCellBlocks();
    std::vector<const GidTypeEntry*> entries(nb);
    std::vector<std::string> mesh_names(nb);
    std::int64_t total_cells = 0;
    for (std::size_t bi = 0; bi < nb; ++bi) {
        const auto cb = rMesh.Cells(bi);
        if (cb.IsRagged() || cb.IsPolyhedron())
            throw WriteError("GiD: cell type '" + cb.Type() +
                             "' has no GiD representation (ragged/polyhedron blocks are "
                             "unsupported)");
        const auto it = gid_type_table().find(cb.Type());
        if (it == gid_type_table().end())
            throw WriteError("GiD: cell type '" + cb.Type() + "' has no verified GiD ordering");
        entries[bi] = &it->second;
        mesh_names[bi] = cb.Type() + "_" + std::to_string(bi);
        total_cells += static_cast<std::int64_t>(cb.NumCells());
    }
    if (static_cast<std::uint64_t>(total_cells) > static_cast<std::uint64_t>(kIntMax))
        throw WriteError("GiD: mesh has too many cells for gidpost's 32-bit API");

    std::vector<std::int64_t> elem_base(nb);
    {
        std::int64_t base = 1;
        for (std::size_t bi = 0; bi < nb; ++bi) {
            elem_base[bi] = base;
            base += static_cast<std::int64_t>(rMesh.Cells(bi).NumCells());
        }
    }

    // A material-id column is written only from one documented key
    // ("gmsh:physical"), only when it is integral and covers every block --
    // never guessed from another convention.
    static const std::string kMatKey = "gmsh:physical";
    const std::string* mat_name = nullptr;
    if (nb > 0 && rMesh.HasCellData(kMatKey) && rMesh.CellDataNumBlocks(kMatKey) == nb) {
        bool all_int = true;
        for (std::size_t bi = 0; bi < nb && all_int; ++bi)
            all_int = gid_is_int_dtype(rMesh.CellData(kMatKey, bi).Dtype());
        if (all_int)
            mat_name = &kMatKey;
    }

    gid_ensure_init();
    const GiD_PostMode gid_mode = gid_post_mode(resolved);

    if (resolved == GidMode::Ascii) {
        const auto paths = gid_ascii_paths(rPath);
        const GiD_FILE fdm = GiD_fOpenPostMeshFile(paths.first.c_str(), gid_mode);
        if (fdm <= 0)
            throw WriteError("GiD: could not open mesh file for writing: " + paths.first);
        gid_write_geometry(fdm, rMesh, entries, mesh_names, mat_name);
        gid_check(GiD_fClosePostMeshFile(fdm), "ClosePostMeshFile");

        const GiD_FILE fdr = GiD_fOpenPostResultFile(paths.second.c_str(), gid_mode);
        if (fdr <= 0)
            throw WriteError("GiD: could not open result file for writing: " + paths.second);
        gid_write_provenance_result(fdr);
        gid_write_point_data(fdr, rMesh, rAnalysisName, stepValue);
        gid_write_cell_data(fdr, rMesh, entries, mesh_names, elem_base, mat_name, rAnalysisName,
                            stepValue);
        gid_check(GiD_fClosePostResultFile(fdr), "ClosePostResultFile");
    } else {
        const GiD_FILE fd = GiD_fOpenPostResultFile(rPath.c_str(), gid_mode);
        if (fd <= 0)
            throw WriteError("GiD: could not open file for writing: " + rPath);
        gid_write_geometry(fd, rMesh, entries, mesh_names, mat_name);
        gid_write_provenance_result(fd);
        gid_write_point_data(fd, rMesh, rAnalysisName, stepValue);
        gid_write_cell_data(fd, rMesh, entries, mesh_names, elem_base, mat_name, rAnalysisName,
                            stepValue);
        gid_check(GiD_fClosePostResultFile(fd), "ClosePostResultFile");
    }
}

}  // namespace meshioplusplus

#else  // !MESHIOPLUSPLUS_HAS_GIDPOST

namespace meshioplusplus {

bool gid_available(GidMode) { return false; }

std::string gid_build_option(GidMode mode) {
    if (mode == GidMode::Hdf5)
        return "-DMESHIOPLUSPLUS_WITH_GIDPOST=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON "
               "-DMESHIOPLUSPLUS_WITH_HDF5=ON";
    return "-DMESHIOPLUSPLUS_WITH_GIDPOST=ON -DMESHIOPLUSPLUS_WITH_ZLIB=ON";
}

// Always defined, so `gid` fails with an error naming the build flags rather
// than a link error or -- worse -- a silent fall-through to another format
// for a `.post.msh` path (the partition_kahip_parts / vtk_codec contract).
void write_gid(const std::string&, const Mesh&, GidMode mode, const std::string&, double) {
    throw WriteError("meshio++: the 'gid' writer needs a build with " + gid_build_option(mode) +
                     " (the vendored gidpost library deflates unconditionally, so zlib is a "
                     "build prerequisite, not an option)");
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_GIDPOST
