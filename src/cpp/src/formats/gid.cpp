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
    // Underscore, not a hyphen: a Julia symbol cannot carry a hyphen, which is
    // why `gradient`'s hyphenated method names needed a translation layer
    // there. Nothing here has to.
    if (rName == "ascii_zipped")
        return GidMode::AsciiZipped;
    throw std::invalid_argument("GiD: unknown mode '" + rName +
                                "' (expected auto, ascii, ascii_zipped, binary, or hdf5)");
}

}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_HAS_GIDPOST

// System includes
#include <cstdint>
#include <limits>
#include <set>
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
// RESOLVED: CIMNE's GiD 6-era hexa20.gif figure numbers the mid-edge nodes
// bottom-ring, VERTICALS, top-ring -- exactly Kratos's INTERNAL order, the
// order Kratos's own GiD writer permutes AWAY from before emitting a file.
// Taken at face value that figure said the identity mapping below is wrong.
// Confirmed that GiD's actual expected order is the one Kratos WRITES, i.e.
// the post-swap order -- which is meshio++'s own hexahedron20 table, hence
// identity here is correct and the figure is outdated. CIMNE's current
// grammar dropped the mid-edge figures entirely for exactly this kind of
// staleness, saying only "hierarchical order ... vertex nodes first, then
// the middle ones" with no order given. The GidOrdering tests below pin that
// no permutation is applied; they do not (and structurally cannot) pin that
// this order is GiD's, which is why this comment states the resolution
// rather than leaving the tests to imply it.
//
// A precision on the evidence, found while deriving hexahedron27/wedge15
// below: the reorder this rests on (`gid_mesh_container.h`) lives ONLY in
// that file's *Conditions*-writing path -- the Elements path writes no
// reorder at all. The conclusion above is unaffected: an Element-agnostic
// Kratos source (`kratos/input_output/vtk_output.cpp`'s Kratos-to-VTK
// conversion, which every Hexahedra3D20 element OR condition goes through
// for VTK/EnSight output, mirrored in `ensight_output.cpp`) independently
// reproduces the identical swap. That broader source is the one this file's
// permutations for hexahedron27/wedge15 actually lean on.
//
// hexahedron27/wedge15 (formerly refused as "not independently verified")
// are now supported via non-identity permutations, `gid_detail::
// gid_cell_perm()` (formats/gid_common.hpp), derived the same way: Kratos's
// own geometry classes (hexahedra_3d_27.h, prism_3d_15.h) order their
// "second ring" of mid-edge nodes as VERTICALS where meshio++'s own table
// (detail/cell_subdivision.cpp) orders it as the opposite tier (top ring /
// top triangle) -- the identical reverse-split pattern hexahedron20 has,
// just with no Conditions-only swap to lean on this time, so these two are
// cross-checked purely against the Element-agnostic vtk_output.cpp/
// ensight_output.cpp source. pyramid13 needs no permutation at all: Kratos's
// own Pyramid3D13 order already matches meshio++'s (vtk_output.cpp
// explicitly skips converting it, i.e. Kratos itself asserts the two agree)
// -- and it is the only one of the three with NO Kratos-GiD precedent
// whatsoever, since Kratos never registers a GiD mesh container for any
// pyramid; its ordering rests on Kratos's internal geometry convention
// alone. See gid_common.hpp's `gid_cell_perm_table()` for the full
// derivation and the permutation arrays themselves.
//
// Anything not in this table -- polygon/polyhedron (GiD has no such type),
// every VTK-Lagrange/higher-degree type -- throws by name rather than
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
        {"hexahedron27", {GiD_Hexahedra, 27}},
        {"wedge", {GiD_Prism, 6}},
        {"wedge15", {GiD_Prism, 15}},
        {"pyramid", {GiD_Pyramid, 5}},
        {"pyramid13", {GiD_Pyramid, 13}},
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
        case GidMode::AsciiZipped:
            return GiD_PostAsciiZipped;
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

        gid_check(
            GiD_fBeginMesh(fd, rMeshNames[bi].c_str(), GiD_3D, entry->mType, static_cast<int>(npc)),
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
        //
        // `perm`, when non-null, is `hexahedron27`/`wedge15`'s node-order
        // permutation (see the cell-type table comment above and
        // gid_common.hpp's derivation): GiD slot j receives meshio++ node
        // perm[j], the exact `dst[c] = src[p[c]]` convention med.cpp's
        // `flatten_f` already uses in this repo.
        const int* perm = gid_detail::gid_cell_perm(cb.Type(), npc);
        std::vector<int> ids(ne);
        std::vector<int> flat_conn(ne * npc);
        for (std::size_t r = 0; r < ne; ++r) {
            ids[r] = static_cast<int>(elem_base + static_cast<std::int64_t>(r));
            for (std::size_t j = 0; j < npc; ++j) {
                const std::size_t src_j = perm ? static_cast<std::size_t>(perm[j]) : j;
                flat_conn[r * npc + j] =
                    static_cast<int>(detail::read_int(conn, r * npc + src_j)) + 1;
            }
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

// meshioplusplus::GidResultType's values are gidpost's GiD_ResultType's, so
// the write below is a static_cast. Pinned here rather than trusted: this is
// the mio_cell_type precedent, and a silent divergence would write a
// plausible file declaring the wrong kind of quantity.
static_assert(static_cast<int>(GidResultType::Scalar) == GiD_Scalar, "GidResultType drift");
static_assert(static_cast<int>(GidResultType::Vector) == GiD_Vector, "GidResultType drift");
static_assert(static_cast<int>(GidResultType::Matrix) == GiD_Matrix, "GidResultType drift");
static_assert(static_cast<int>(GidResultType::PlainDeformationMatrix) == GiD_PlainDeformationMatrix,
              "GidResultType drift");
static_assert(static_cast<int>(GidResultType::MainMatrix) == GiD_MainMatrix, "GidResultType drift");
static_assert(static_cast<int>(GidResultType::LocalAxes) == GiD_LocalAxes, "GidResultType drift");
static_assert(static_cast<int>(GidResultType::ComplexScalar) == GiD_ComplexScalar,
              "GidResultType drift");
static_assert(static_cast<int>(GidResultType::ComplexVector) == GiD_ComplexVector,
              "GidResultType drift");
static_assert(static_cast<int>(GidResultType::ComplexMatrix) == GiD_ComplexMatrix,
              "GidResultType drift");

/**
 * @brief The `GidResultType` a caller declared for @p rName, if any.
 *
 * The declaration rides a `field_data` entry (`kGidResultTypePrefix + name`)
 * rather than a side-channel struct so that it reaches every surface: the
 * registry's writers take `(path, mesh)` and could not carry a struct, which
 * would have left the CLI, WASM, the C API and Fortran unable to declare
 * anything. See `gid.hpp`.
 *
 * @throws WriteError when the value is out of range, or when @p k is not a
 *         count that type accepts -- never a silent fallback to splitting,
 *         which would quietly ignore an explicit request (write_options.hpp's
 *         standing rule).
 */
const gid_detail::GidResultTypeEntry* gid_declared_result_type(const Mesh& rMesh,
                                                               const std::string& rName,
                                                               std::size_t k) {
    const std::string key = std::string(kGidResultTypePrefix) + rName;
    if (!rMesh.HasFieldData(key))
        return nullptr;

    const NDArray& decl = rMesh.FieldData(key);
    if (decl.Size() == 0)
        throw WriteError("GiD: field_data['" + key + "'] is empty; it must hold one integer " +
                         "GidResultType value");
    const auto raw = static_cast<std::int64_t>(detail::read_double(decl, 0));
    const gid_detail::GidResultTypeEntry* entry =
        gid_detail::gid_find_result_type(static_cast<GidResultType>(raw));
    if (entry == nullptr)
        throw WriteError("GiD: field_data['" + key + "'] declares result type " +
                         std::to_string(raw) + ", which is not a GidResultType value (0..8)");
    if (!gid_result_dim_is_legal(entry->mType, k))
        throw WriteError("GiD: '" + rName + "' is declared " + entry->mName + " but has " +
                         std::to_string(k) + " components; " + entry->mName + " accepts " +
                         gid_detail::gid_legal_dims_text(*entry));
    return entry;
}

/**
 * @brief The Gauss-point count declared for @p rName, or 1.
 *
 * `field_data[kGidGaussPointsPrefix + name]` -- see gid.hpp. Absent means the
 * historical one-value-per-element behaviour, which is why an ordinary mesh
 * never touches any of the G>1 machinery below.
 *
 * @throws WriteError on a non-positive G, or a G that does not divide the
 *         array's own column count (which would mean the flat `(ncells, G*k)`
 *         layout cannot be split into whole components).
 */
std::size_t gid_declared_gauss_points(const Mesh& rMesh, const std::string& rName,
                                      std::size_t cols) {
    const std::string key = std::string(kGidGaussPointsPrefix) + rName;
    if (!rMesh.HasFieldData(key))
        return 1;

    const NDArray& decl = rMesh.FieldData(key);
    if (decl.Size() == 0)
        throw WriteError("GiD: field_data['" + key +
                         "'] is empty; it must hold one integer Gauss-point count");
    const auto raw = static_cast<std::int64_t>(detail::read_double(decl, 0));
    if (raw < 1)
        throw WriteError("GiD: field_data['" + key + "'] declares " + std::to_string(raw) +
                         " Gauss points; the count must be at least 1");
    const auto g = static_cast<std::size_t>(raw);
    if (cols % g != 0)
        throw WriteError("GiD: '" + rName + "' is declared " + std::to_string(g) +
                         " Gauss points but has " + std::to_string(cols) + " columns, which " +
                         std::to_string(g) +
                         " does not divide; the layout is (ncells, G*components)");
    return g;
}

/**
 * @brief Declares one Gauss-point set, `Internal` or `Given`.
 *
 * `Internal` (GiD places the points) is only legal for specific counts per
 * element family; any other G must list its natural coordinates explicitly.
 * Supplying coordinates for a count GiD *could* have placed is honoured --
 * a solver's own quadrature rule need not match GiD's -- but omitting them
 * for one it cannot is an error naming the legal counts rather than a file
 * GiD would silently reject.
 */
void gid_declare_gauss_set(GiD_FILE fd, const Mesh& rMesh, const std::string& rGaussName,
                           const std::string& rMeshioType, GiD_ElementType EType,
                           const std::string& rMeshName, std::size_t g) {
    const gid_detail::GidGaussCountEntry* fam = gid_detail::gid_find_gauss_counts(rMeshioType);
    const std::string coords_key = gid_detail::gid_gauss_coords_key(rMeshioType, g);
    const bool has_coords = rMesh.HasFieldData(coords_key);

    if (has_coords && fam != nullptr && fam->mCoordDim == 0)
        throw WriteError("GiD: '" + rMeshioType +
                         "' cannot use given Gauss-point coordinates (GiD forbids "
                         "'Natural Coordinates: Given' for line elements); remove field_data['" +
                         coords_key + "']");
    if (!has_coords && !gid_detail::gid_gauss_count_is_internal(rMeshioType, g))
        throw WriteError("GiD: '" + rMeshioType + "' cannot place " + std::to_string(g) +
                         " Gauss points itself (it accepts " +
                         gid_detail::gid_internal_counts_text(*fam) +
                         "); supply their natural coordinates in field_data['" + coords_key + "']");

    const std::size_t dim = fam != nullptr ? fam->mCoordDim : 3;
    if (has_coords) {
        const NDArray& coords = rMesh.FieldData(coords_key);
        if (coords.Size() != g * dim)
            throw WriteError("GiD: field_data['" + coords_key + "'] has " +
                             std::to_string(coords.Size()) + " values; " + std::to_string(g) +
                             " points on '" + rMeshioType + "' need " + std::to_string(g * dim) +
                             " (" + std::to_string(dim) + " per point)");
    }

    gid_check(GiD_fBeginGaussPoint(fd, rGaussName.c_str(), EType, rMeshName.c_str(),
                                   static_cast<int>(g), /*NodesIncluded=*/0,
                                   /*InternalCoord=*/has_coords ? 0 : 1),
              "BeginGaussPoint('" + rGaussName + "')");
    if (has_coords) {
        const NDArray& coords = rMesh.FieldData(coords_key);
        for (std::size_t i = 0; i < g; ++i) {
            const double a = detail::read_double(coords, i * dim + 0);
            const double b = detail::read_double(coords, i * dim + 1);
            if (dim == 2)
                gid_check(GiD_fWriteGaussPoint2D(fd, a, b), "WriteGaussPoint2D");
            else
                gid_check(
                    GiD_fWriteGaussPoint3D(fd, a, b, detail::read_double(coords, i * dim + 2)),
                    "WriteGaussPoint3D");
        }
    }
    gid_check(GiD_fEndGaussPoint(fd), "EndGaussPoint('" + rGaussName + "')");
}

/// The Gauss-point set name for a block at @p g points. `G == 1` keeps the
/// historical `gp_<mesh>` spelling EXACTLY -- byte-identity for every mesh
/// that never asks for multiple Gauss points depends on it.
std::string gid_gauss_set_name(const std::string& rMeshName, std::size_t g) {
    return g == 1 ? "gp_" + rMeshName : "gp" + std::to_string(g) + "_" + rMeshName;
}

// GiD_ResultType's valid component counts are irregular (Scalar 1; Vector
// 2/3/4; Matrix 3/6; MainMatrix 12; ...) and gidpost does not validate an
// unsupported count itself -- it silently emits a malformed file. A declared
// type is validated above and written as-is; with no declaration only 1
// (Scalar) and 2/3 (Vector) are mapped directly and anything else splits into
// that many named scalars, because a bare component count is genuinely
// ambiguous (see gid.hpp's k-component rule).
void gid_write_result_array(GiD_FILE fd, const Mesh& rMesh, const NDArray& rArr, std::size_t rows,
                            const std::vector<int>& rIds, const std::string& rName,
                            const std::string& rAnalysis, double step, GiD_ResultLocation loc,
                            const std::string& rGaussName, std::size_t g = 1) {
    // The array is (rows, G*k) laid out Gauss-point-major, so the COMPONENT
    // count -- what a GiD result type's legal widths are checked against, and
    // what goes in the Values row -- is cols/G, not cols. Passing cols here
    // would misdeclare every G>1 array's type and width.
    const std::size_t cols = rArr.Shape().size() >= 2 ? rArr.Shape()[1] : 1;
    const std::size_t k = cols / g;
    const char* gauss = rGaussName.empty() ? nullptr : rGaussName.c_str();
    const gid_detail::GidResultTypeEntry* declared = gid_declared_result_type(rMesh, rName, k);

    // GiD wants G rows per element; gidpost's own writer suppresses a repeated
    // id, so handing it each id G times produces exactly the compact form the
    // reader expects, with no id bookkeeping here.
    const std::size_t out_rows = rows * g;
    std::vector<int> ids_g;
    const std::vector<int>* ids = &rIds;
    if (g != 1) {
        ids_g.reserve(out_rows);
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t p = 0; p < g; ++p)
                ids_g.push_back(rIds[r]);
        ids = &ids_g;
    }

    // ONE inference rule, shared with the reader (which needs it to decide
    // whether a file's declared type carries information worth recording).
    // Two expressions of it would drift, and the drift is silent: the writer
    // would emit a type the reader then declines to record, so the round trip
    // would quietly lose the declaration.
    const gid_detail::GidResultTypeEntry* chosen =
        declared != nullptr ? declared : gid_detail::gid_inferred_result_type(k);

    if (chosen != nullptr) {
        const auto rtype = static_cast<GiD_ResultType>(chosen->mType);
        // Source and destination share the Gauss-point-major layout, so this
        // is a straight copy for any G -- the (r, gp, c) index is r*cols +
        // p*k + c on both sides.
        std::vector<double> vals(out_rows * k);
        for (std::size_t i = 0; i < out_rows * k; ++i)
            vals[i] = detail::read_double(rArr, i);
        gid_check(
            GiD_fWriteResultBlock(fd, rName.c_str(), rAnalysis.c_str(), step, rtype, loc, gauss,
                                  nullptr, 0, nullptr, nullptr, static_cast<int>(out_rows),
                                  ids->data(), static_cast<int>(k), vals.data()),
            "WriteResultBlock('" + rName + "')");
        return;
    }

    detail::provenance_note("result-split", rName + ": k=" + std::to_string(k) +
                                                " has no GiD result type; written as " +
                                                std::to_string(k) + " scalars");
    std::vector<double> col(out_rows);
    for (std::size_t c = 0; c < k; ++c) {
        for (std::size_t i = 0; i < out_rows; ++i)
            col[i] = detail::read_double(rArr, i * k + c);
        const std::string cname = rName + "_" + std::to_string(c + 1);
        gid_check(GiD_fWriteResultBlock(fd, cname.c_str(), rAnalysis.c_str(), step, GiD_Scalar, loc,
                                        gauss, nullptr, 0, nullptr, nullptr,
                                        static_cast<int>(out_rows), ids->data(), 1, col.data()),
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
        gid_write_result_array(fd, rMesh, rMesh.PointData(name), np, ids, name, rAnalysis, step,
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

    // Gauss-point sets are keyed by (block, G), not by block: a set's identity
    // in GiD IS its point count, so two arrays on one block that declare
    // different counts need two sets, while two declaring the same count share
    // one. Collected first so each set is declared exactly once, before any
    // result references it.
    //
    // Every block still gets its G=1 set unconditionally, under the historical
    // "gp_<mesh>" name and in the historical order, so a mesh with no
    // Gauss-point declarations emits byte-identical bytes to before this
    // existed.
    std::vector<std::set<std::size_t>> block_counts(nb);
    for (std::size_t bi = 0; bi < nb; ++bi)
        block_counts[bi].insert(1);
    for (const auto& name : rMesh.CellDataNames()) {
        if (pMatName != nullptr && name == *pMatName)
            continue;
        if (rMesh.CellDataNumBlocks(name) != nb)
            continue;  // reported below, where the array is written
        for (std::size_t bi = 0; bi < nb; ++bi) {
            const NDArray& arr = rMesh.CellData(name, bi);
            const std::size_t cols = arr.Shape().size() >= 2 ? arr.Shape()[1] : 1;
            block_counts[bi].insert(gid_declared_gauss_points(rMesh, name, cols));
        }
    }

    for (std::size_t bi = 0; bi < nb; ++bi)
        for (std::size_t g : block_counts[bi])
            gid_declare_gauss_set(fd, rMesh, gid_gauss_set_name(rMeshNames[bi], g),
                                  rMesh.Cells(bi).Type(), rEntries[bi]->mType, rMeshNames[bi], g);

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
            const NDArray& arr = rMesh.CellData(name, bi);
            const std::size_t cols = arr.Shape().size() >= 2 ? arr.Shape()[1] : 1;
            const std::size_t g = gid_declared_gauss_points(rMesh, name, cols);
            gid_write_result_array(fd, rMesh, arr, ne, ids, name, rAnalysis, step,
                                   GiD_OnGaussPoints, gid_gauss_set_name(rMeshNames[bi], g), g);
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

    // AsciiZipped is the same two-sibling-file layout as Ascii -- it is the
    // identical text, only gzipped -- so it shares this branch entirely.
    if (resolved == GidMode::Ascii || resolved == GidMode::AsciiZipped) {
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

bool gid_available(GidMode) {
    return false;
}

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
