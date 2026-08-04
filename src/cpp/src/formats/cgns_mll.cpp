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
// The optional cgnslib (CGNS Mid-Level Library) reader. Additive: the
// hand-rolled ADF-over-HDF5 reader in cgns.cpp is untouched and remains the
// default, so a build without this flag behaves exactly as before. See
// formats/cgns.hpp's read_cgns_mll doc comment for what it buys.
//
// The whole file is behind the guard, so CMake's GLOB_RECURSE still compiles it
// down to an empty translation unit when the flag is off -- the same shape
// cgns.cpp uses for MESHIOPLUSPLUS_HAS_HDF5.

// Project includes
#include "meshioplusplus/formats/cgns.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/format_compat.hpp"

#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB
#include <cgnslib.h>
#endif

namespace meshioplusplus {

bool cgns_has_cgnslib() {
#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB
    return true;
#else
    return false;
#endif
}

#ifndef MESHIOPLUSPLUS_HAS_CGNSLIB

Mesh read_cgns_mll(const std::string& rPath) {
    // Always present and throwing by name -- the partition_kahip_parts
    // contract. A link error would break the Python-fallback contract, and a
    // silent downgrade to the raw-HDF5 reader would answer a question the
    // caller did not ask (that reader cannot open an ADF file at all).
    throw ReadError(detail::format_compat(
        "meshio++: cannot read '{}' through cgnslib: this build has no cgnslib support "
        "(rebuild with -DMESHIOPLUSPLUS_WITH_CGNSLIB=ON and CGNS_ROOT pointing at an install)",
        rPath));
}

#else

namespace {

/// Raise with cgnslib's own error text, which is the only useful diagnostic it
/// gives -- the return code alone says nothing about what was wrong.
[[noreturn]] void cgns_mll_fail(const std::string& rWhat, const std::string& rPath) {
    throw ReadError(detail::format_compat("meshio++: CGNS (cgnslib): {} while reading '{}': {}",
                                          rWhat, rPath, cg_get_error()));
}

/// RAII for the library's integer file handle: every early return below would
/// otherwise leak it, and cgnslib keeps global state per open file.
class CgnsFile {
public:
    CgnsFile(const std::string& rPath) : mPath(rPath) {
        if (cg_open(rPath.c_str(), CG_MODE_READ, &mFn) != CG_OK)
            cgns_mll_fail("cg_open failed", rPath);
    }
    ~CgnsFile() {
        if (mFn >= 0)
            cg_close(mFn);
    }
    CgnsFile(const CgnsFile&) = delete;
    CgnsFile& operator=(const CgnsFile&) = delete;
    int Fn() const { return mFn; }

private:
    std::string mPath;
    int mFn = -1;
};

/// meshio++ name for a fixed-size CGNS ElementType_t, or empty when the type is
/// one this reader deliberately does not claim (see the header).
std::string cgns_mll_meshio_name(CGNS_ENUMT(ElementType_t) type) {
    switch (type) {
        case CGNS_ENUMV(NODE):
            return "vertex";
        case CGNS_ENUMV(BAR_2):
            return "line";
        case CGNS_ENUMV(BAR_3):
            return "line3";
        case CGNS_ENUMV(TRI_3):
            return "triangle";
        case CGNS_ENUMV(TRI_6):
            return "triangle6";
        case CGNS_ENUMV(QUAD_4):
            return "quad";
        case CGNS_ENUMV(QUAD_8):
            return "quad8";
        case CGNS_ENUMV(QUAD_9):
            return "quad9";
        case CGNS_ENUMV(TETRA_4):
            return "tetra";
        case CGNS_ENUMV(TETRA_10):
            return "tetra10";
        case CGNS_ENUMV(PYRA_5):
            return "pyramid";
        case CGNS_ENUMV(PYRA_13):
            return "pyramid13";
        case CGNS_ENUMV(PYRA_14):
            return "pyramid14";
        case CGNS_ENUMV(PENTA_6):
            return "wedge";
        case CGNS_ENUMV(PENTA_15):
            return "wedge15";
        case CGNS_ENUMV(PENTA_18):
            return "wedge18";
        case CGNS_ENUMV(HEXA_8):
            return "hexahedron";
        case CGNS_ENUMV(HEXA_20):
            return "hexahedron20";
        case CGNS_ENUMV(HEXA_27):
            return "hexahedron27";
        default:
            return std::string();
    }
}

/// The SIDS<->meshio node permutation for the types whose orderings differ.
/// Self-inverse, so one table serves both directions -- the same tables
/// cgns.cpp's cgns_type_table() carries, restated here rather than exported
/// because that one is file-private and its shape (name + code + perm) does not
/// fit a lookup keyed on the MLL's enum.
const std::vector<int>* cgns_mll_perm(const std::string& rName) {
    static const std::map<std::string, std::vector<int> > perms = {
        {"wedge15", {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}},
        {"wedge18", {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17}},
        {"hexahedron20", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}},
        {"hexahedron27", {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 16, 17,
                          18, 19, 12, 13, 14, 15, 24, 22, 21, 23, 20, 25, 26}},
    };
    auto it = perms.find(rName);
    return it == perms.end() ? nullptr : &it->second;
}

/// One section's raw description, gathered before anything is decoded.
struct MllSection {
    int mIndex = 0;
    std::string mName;
    CGNS_ENUMT(ElementType_t) mType = CGNS_ENUMV(ElementTypeNull);
    cgsize_t mStart = 0;
    cgsize_t mEnd = 0;
};

/// Read one fixed-size section into a rectangular block.
void cgns_mll_read_fixed(int fn, int B, int Z, const MllSection& rSec, const std::string& rMeshio,
                         const std::string& rPath, Mesh& rMesh) {
    const int npc = [&] {
        int n = 0;
        cg_npe(rSec.mType, &n);
        return n;
    }();
    if (npc <= 0)
        cgns_mll_fail("cg_npe returned no node count for section '" + rSec.mName + "'", rPath);

    const std::size_t ncells = static_cast<std::size_t>(rSec.mEnd - rSec.mStart + 1);
    std::vector<cgsize_t> raw(ncells * static_cast<std::size_t>(npc));
    if (cg_elements_read(fn, B, Z, rSec.mIndex, raw.data(), nullptr) != CG_OK)
        cgns_mll_fail("cg_elements_read failed for section '" + rSec.mName + "'", rPath);

    NDArray conn = NDArray::Uninit(DType::Int64, {ncells, static_cast<std::size_t>(npc)});
    std::int64_t* dst = conn.As<std::int64_t>();
    const std::vector<int>* perm = cgns_mll_perm(rMeshio);
    for (std::size_t c = 0; c < ncells; ++c) {
        for (int k = 0; k < npc; ++k) {
            // CGNS node ids are 1-based; the permutation is self-inverse, so
            // the same scatter serves read and write.
            const std::int64_t v =
                static_cast<std::int64_t>(raw[c * static_cast<std::size_t>(npc) + k]) - 1;
            dst[c * static_cast<std::size_t>(npc) + (perm ? (*perm)[k] : k)] = v;
        }
    }
    rMesh.AddCellBlock(rMeshio, std::move(conn));
}

/// A face list decoded from an NGON_n section: CSR over 0-based node ids.
struct MllFaces {
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mOffsets;  // numFaces + 1
    cgsize_t mFirstId = 0;               // the section's ElementRange start
    std::size_t Count() const { return mOffsets.empty() ? 0 : mOffsets.size() - 1; }
};

MllFaces cgns_mll_read_ngon(int fn, int B, int Z, const MllSection& rSec,
                            const std::string& rPath) {
    const std::size_t nfaces = static_cast<std::size_t>(rSec.mEnd - rSec.mStart + 1);
    cgsize_t data_size = 0;
    if (cg_ElementDataSize(fn, B, Z, rSec.mIndex, &data_size) != CG_OK)
        cgns_mll_fail("cg_ElementDataSize failed for NGON_n section '" + rSec.mName + "'", rPath);

    std::vector<cgsize_t> elems(static_cast<std::size_t>(data_size));
    std::vector<cgsize_t> offsets(nfaces + 1, 0);
    // cg_poly_elements_read is what makes this worth doing at all: it presents
    // both the CGNS 3.x layout (each face prefixed by its node count) and the
    // 4.0 one (a separate ElementStartOffset array) through this one shape.
    if (cg_poly_elements_read(fn, B, Z, rSec.mIndex, elems.data(), offsets.data(), nullptr) !=
        CG_OK)
        cgns_mll_fail("cg_poly_elements_read failed for NGON_n section '" + rSec.mName + "'",
                      rPath);

    MllFaces out;
    out.mFirstId = rSec.mStart;
    out.mOffsets.reserve(nfaces + 1);
    out.mNodes.reserve(static_cast<std::size_t>(data_size));
    out.mOffsets.push_back(0);
    for (std::size_t f = 0; f < nfaces; ++f) {
        for (cgsize_t i = offsets[f]; i < offsets[f + 1]; ++i)
            out.mNodes.push_back(static_cast<std::int64_t>(elems[static_cast<std::size_t>(i)]) - 1);
        out.mOffsets.push_back(static_cast<std::int64_t>(out.mNodes.size()));
    }
    return out;
}

/// Rejoin `<base>_0.._k-1` siblings into one k-component array, mirroring
/// cgns.cpp's cgns_read_solution. Names with no such contiguous run stay scalar.
void cgns_mll_group_components(const std::vector<std::string>& rNames,
                               std::vector<std::string>& rBases,
                               std::vector<std::vector<std::size_t> >& rGroups) {
    std::vector<char> taken(rNames.size(), 0);
    for (std::size_t i = 0; i < rNames.size(); ++i) {
        if (taken[i])
            continue;
        const std::string& n = rNames[i];
        const std::size_t us = n.rfind('_');
        bool grouped = false;
        if (us != std::string::npos && us + 1 < n.size() &&
            n.find_first_not_of("0123456789", us + 1) == std::string::npos &&
            n.substr(us + 1) == "0") {
            const std::string base = n.substr(0, us);
            std::vector<std::size_t> run{i};
            for (std::size_t k = 1;; ++k) {
                const std::string want = base + "_" + std::to_string(k);
                std::size_t found = rNames.size();
                for (std::size_t j = 0; j < rNames.size(); ++j)
                    if (!taken[j] && j != i && rNames[j] == want)
                        found = j;
                if (found == rNames.size())
                    break;
                run.push_back(found);
            }
            if (run.size() > 1) {
                for (std::size_t j : run)
                    taken[j] = 1;
                rBases.push_back(base);
                rGroups.push_back(std::move(run));
                grouped = true;
            }
        }
        if (!grouped) {
            taken[i] = 1;
            rBases.push_back(n);
            rGroups.push_back({i});
        }
    }
}

/// Read every FlowSolution_t into point_data / cell_data.
void cgns_mll_read_solutions(int fn, int B, int Z, std::size_t NumPoints, Mesh& rMesh,
                             const std::string& rPath) {
    int nsols = 0;
    if (cg_nsols(fn, B, Z, &nsols) != CG_OK || nsols < 1)
        return;

    // cell_data is per-BLOCK in meshio++ but per-ZONE in CGNS, so a zone-wide
    // array is split back across the blocks in ElementRange order -- the same
    // rule cgns.cpp applies, and the reason a mixed-dimension mesh cannot carry
    // it at all.
    std::vector<std::size_t> block_cells;
    for (const auto cb : rMesh.CellRange())
        block_cells.push_back(cb.NumCells());
    std::size_t total_cells = 0;
    for (std::size_t n : block_cells)
        total_cells += n;

    for (int S = 1; S <= nsols; ++S) {
        char sol_name[33] = {0};
        CGNS_ENUMT(GridLocation_t) loc = CGNS_ENUMV(GridLocationNull);
        if (cg_sol_info(fn, B, Z, S, sol_name, &loc) != CG_OK)
            cgns_mll_fail("cg_sol_info failed", rPath);
        if (loc != CGNS_ENUMV(Vertex) && loc != CGNS_ENUMV(CellCenter)) {
            log::warn(
                "CGNS (cgnslib): FlowSolution '{}' has GridLocation {} (only Vertex and "
                "CellCenter are mapped); skipping it",
                sol_name, static_cast<int>(loc));
            continue;
        }
        const bool vertex = loc == CGNS_ENUMV(Vertex);
        const std::size_t rows = vertex ? NumPoints : total_cells;
        if (rows == 0)
            continue;

        int nfields = 0;
        if (cg_nfields(fn, B, Z, S, &nfields) != CG_OK)
            cgns_mll_fail("cg_nfields failed", rPath);
        std::vector<std::string> names;
        for (int F = 1; F <= nfields; ++F) {
            CGNS_ENUMT(DataType_t) dt = CGNS_ENUMV(DataTypeNull);
            char fname[33] = {0};
            if (cg_field_info(fn, B, Z, S, F, &dt, fname) != CG_OK)
                cgns_mll_fail("cg_field_info failed", rPath);
            names.emplace_back(fname);
        }

        std::vector<std::string> bases;
        std::vector<std::vector<std::size_t> > groups;
        cgns_mll_group_components(names, bases, groups);

        const cgsize_t rmin = 1;
        const cgsize_t rmax = static_cast<cgsize_t>(rows);
        for (std::size_t g = 0; g < groups.size(); ++g) {
            const std::size_t ncomp = groups[g].size();
            std::vector<double> buf(rows);
            NDArray arr = ncomp == 1 ? NDArray::Uninit(DType::Float64, {rows})
                                     : NDArray::Uninit(DType::Float64, {rows, ncomp});
            double* dst = arr.As<double>();
            for (std::size_t k = 0; k < ncomp; ++k) {
                if (cg_field_read(fn, B, Z, S, names[groups[g][k]].c_str(), CGNS_ENUMV(RealDouble),
                                  &rmin, &rmax, buf.data()) != CG_OK)
                    cgns_mll_fail("cg_field_read failed for '" + names[groups[g][k]] + "'", rPath);
                for (std::size_t r = 0; r < rows; ++r)
                    dst[r * ncomp + k] = buf[r];
            }
            if (vertex) {
                rMesh.AddPointData(bases[g], std::move(arr));
            } else {
                std::vector<NDArray> blocks;
                std::size_t at = 0;
                for (std::size_t n : block_cells) {
                    NDArray b = ncomp == 1 ? NDArray::Uninit(DType::Float64, {n})
                                           : NDArray::Uninit(DType::Float64, {n, ncomp});
                    std::memcpy(b.Data(), arr.As<double>() + at * ncomp,
                                n * ncomp * sizeof(double));
                    blocks.push_back(std::move(b));
                    at += n;
                }
                rMesh.AddCellData(bases[g], std::move(blocks));
            }
        }
    }
}

}  // namespace

Mesh read_cgns_mll(const std::string& rPath) {
    CgnsFile file(rPath);
    const int fn = file.Fn();

    int nbases = 0;
    if (cg_nbases(fn, &nbases) != CG_OK || nbases < 1)
        cgns_mll_fail("no CGNSBase_t found", rPath);
    if (nbases > 1)
        log::warn("CGNS (cgnslib): '{}' has {} bases; reading the first only", rPath, nbases);
    const int B = 1;

    char base_name[33] = {0};
    int cell_dim = 0, phys_dim = 0;
    if (cg_base_read(fn, B, base_name, &cell_dim, &phys_dim) != CG_OK)
        cgns_mll_fail("cg_base_read failed", rPath);

    int nzones = 0;
    if (cg_nzones(fn, B, &nzones) != CG_OK || nzones < 1)
        cgns_mll_fail("no Zone_t found", rPath);
    if (nzones > 1)
        log::warn("CGNS (cgnslib): '{}' has {} zones; reading the first only", rPath, nzones);
    const int Z = 1;

    CGNS_ENUMT(ZoneType_t) ztype = CGNS_ENUMV(ZoneTypeNull);
    if (cg_zone_type(fn, B, Z, &ztype) != CG_OK)
        cgns_mll_fail("cg_zone_type failed", rPath);
    if (ztype != CGNS_ENUMV(Unstructured))
        throw ReadError(detail::format_compat(
            "meshio++: CGNS (cgnslib): '{}' zone 1 is Structured; meshio++ reads unstructured "
            "zones only",
            rPath));

    char zone_name[33] = {0};
    cgsize_t zsize[9] = {0};
    if (cg_zone_read(fn, B, Z, zone_name, zsize) != CG_OK)
        cgns_mll_fail("cg_zone_read failed", rPath);
    const std::size_t npoints = static_cast<std::size_t>(zsize[0]);

    // --- coordinates -------------------------------------------------------
    int ncoords = 0;
    if (cg_ncoords(fn, B, Z, &ncoords) != CG_OK)
        cgns_mll_fail("cg_ncoords failed", rPath);
    const std::size_t dim = static_cast<std::size_t>(std::max(2, std::min(3, ncoords)));

    Mesh mesh;
    {
        NDArray points = NDArray::Uninit(DType::Float64, {npoints, dim});
        double* dst = points.As<double>();
        std::vector<double> buf(npoints);
        const cgsize_t rmin = 1;
        const cgsize_t rmax = static_cast<cgsize_t>(npoints);
        static const char* kNames[3] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
        for (std::size_t d = 0; d < dim; ++d) {
            if (cg_coord_read(fn, B, Z, kNames[d], CGNS_ENUMV(RealDouble), &rmin, &rmax,
                              buf.data()) != CG_OK)
                cgns_mll_fail(std::string("cg_coord_read failed for ") + kNames[d], rPath);
            for (std::size_t i = 0; i < npoints; ++i)
                dst[i * dim + d] = buf[i];
        }
        mesh.AssignPoints(std::move(points));
    }

    // --- sections ----------------------------------------------------------
    int nsections = 0;
    if (cg_nsections(fn, B, Z, &nsections) != CG_OK)
        cgns_mll_fail("cg_nsections failed", rPath);

    std::vector<MllSection> secs;
    for (int S = 1; S <= nsections; ++S) {
        MllSection sec;
        sec.mIndex = S;
        char name[33] = {0};
        int nbndry = 0, parent_flag = 0;
        if (cg_section_read(fn, B, Z, S, name, &sec.mType, &sec.mStart, &sec.mEnd, &nbndry,
                            &parent_flag) != CG_OK)
            cgns_mll_fail("cg_section_read failed", rPath);
        sec.mName = name;
        secs.push_back(std::move(sec));
    }
    // Ascending ElementRange, which reproduces a writer's own block order.
    std::sort(secs.begin(), secs.end(),
              [](const MllSection& a, const MllSection& b) { return a.mStart < b.mStart; });

    // NGON_n/NFACE_n are read together: NFACE_n's cells reference face ids that
    // only NGON_n can resolve, so both are gathered before either is emitted.
    std::vector<MllFaces> ngons;
    std::vector<const MllSection*> nfaces;
    for (const MllSection& sec : secs) {
        if (sec.mType == CGNS_ENUMV(NGON_n)) {
            ngons.push_back(cgns_mll_read_ngon(fn, B, Z, sec, rPath));
            continue;
        }
        if (sec.mType == CGNS_ENUMV(NFACE_n)) {
            nfaces.push_back(&sec);
            continue;
        }
        if (sec.mType == CGNS_ENUMV(MIXED))
            throw ReadError(detail::format_compat(
                "meshio++: CGNS (cgnslib): element section '{}' has ElementType MIXED (20); "
                "MIXED sections are not supported.",
                sec.mName));
        const std::string meshio = cgns_mll_meshio_name(sec.mType);
        if (meshio.empty())
            throw ReadError(detail::format_compat(
                "meshio++: CGNS (cgnslib): element section '{}' has ElementType {}, whose node "
                "ordering meshio++ has not verified (the cubic/quartic Lagrange family); "
                "refusing to guess.",
                sec.mName, static_cast<int>(sec.mType)));
        cgns_mll_read_fixed(fn, B, Z, sec, meshio, rPath, mesh);
    }

    if (!ngons.empty() && nfaces.empty()) {
        // Faces with no cells referencing them: a face mesh, not a volume one.
        for (const MllFaces& faces : ngons) {
            std::vector<std::vector<std::int64_t> > rows(faces.Count());
            for (std::size_t f = 0; f < faces.Count(); ++f)
                rows[f].assign(faces.mNodes.begin() + faces.mOffsets[f],
                               faces.mNodes.begin() + faces.mOffsets[f + 1]);
            if (!rows.empty())
                mesh.AddPolygonBlock("polygon", std::move(rows));
        }
    } else if (!nfaces.empty()) {
        if (ngons.empty())
            throw ReadError(detail::format_compat(
                "meshio++: CGNS (cgnslib): '{}' has an NFACE_n section but no NGON_n to resolve "
                "its face ids against",
                rPath));
        // One flat face table across every NGON_n section, indexed by the CGNS
        // element id the NFACE_n entries carry.
        std::map<cgsize_t, std::pair<const MllFaces*, std::size_t> > by_id;
        for (const MllFaces& faces : ngons)
            for (std::size_t f = 0; f < faces.Count(); ++f)
                by_id.emplace(faces.mFirstId + static_cast<cgsize_t>(f), std::make_pair(&faces, f));

        for (const MllSection* sec : nfaces) {
            const std::size_t ncells = static_cast<std::size_t>(sec->mEnd - sec->mStart + 1);
            cgsize_t data_size = 0;
            if (cg_ElementDataSize(fn, B, Z, sec->mIndex, &data_size) != CG_OK)
                cgns_mll_fail("cg_ElementDataSize failed for NFACE_n '" + sec->mName + "'", rPath);
            std::vector<cgsize_t> elems(static_cast<std::size_t>(data_size));
            std::vector<cgsize_t> offsets(ncells + 1, 0);
            if (cg_poly_elements_read(fn, B, Z, sec->mIndex, elems.data(), offsets.data(),
                                      nullptr) != CG_OK)
                cgns_mll_fail("cg_poly_elements_read failed for NFACE_n '" + sec->mName + "'",
                              rPath);

            // Group by unique node count into polyhedron<N>, the convention the
            // OpenFOAM and EnSight readers already use.
            std::vector<std::vector<std::vector<std::int64_t> > > cells(ncells);
            std::vector<std::size_t> node_counts(ncells, 0);
            for (std::size_t c = 0; c < ncells; ++c) {
                std::vector<std::int64_t> uniq;
                for (cgsize_t i = offsets[c]; i < offsets[c + 1]; ++i) {
                    const cgsize_t signed_id = elems[static_cast<std::size_t>(i)];
                    const cgsize_t id = signed_id < 0 ? -signed_id : signed_id;
                    auto it = by_id.find(id);
                    if (it == by_id.end())
                        throw ReadError(detail::format_compat(
                            "meshio++: CGNS (cgnslib): NFACE_n section '{}' references face id "
                            "{}, which no NGON_n section defines",
                            sec->mName, static_cast<long long>(id)));
                    const MllFaces& faces = *it->second.first;
                    const std::size_t f = it->second.second;
                    std::vector<std::int64_t> ring(faces.mNodes.begin() + faces.mOffsets[f],
                                                   faces.mNodes.begin() + faces.mOffsets[f + 1]);
                    // A negative id means "this face, traversed the other way"
                    // -- CGNS's way of orienting a shared face outward from
                    // each of the two cells that use it.
                    if (signed_id < 0)
                        std::reverse(ring.begin(), ring.end());
                    uniq.insert(uniq.end(), ring.begin(), ring.end());
                    cells[c].push_back(std::move(ring));
                }
                std::sort(uniq.begin(), uniq.end());
                uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                node_counts[c] = uniq.size();
            }

            std::vector<std::size_t> order;
            std::map<std::size_t, std::vector<std::size_t> > groups;
            for (std::size_t c = 0; c < ncells; ++c) {
                if (groups.find(node_counts[c]) == groups.end())
                    order.push_back(node_counts[c]);
                groups[node_counts[c]].push_back(c);
            }
            for (std::size_t n : order) {
                std::vector<std::vector<std::vector<std::int64_t> > > group;
                group.reserve(groups[n].size());
                for (std::size_t c : groups[n])
                    group.push_back(std::move(cells[c]));
                mesh.AddPolyhedronBlock("polyhedron" + std::to_string(n), std::move(group));
            }
        }
    }

    // --- FlowSolution_t ----------------------------------------------------
    //
    // CGNS has no component concept, so a k-component meshio++ array is stored
    // as k siblings suffixed `_0.._k-1` and rejoined here from a CONTIGUOUS
    // run starting at 0 -- exactly the convention cgns.cpp writes and
    // documents. Anything else (a lone `foo_7`, a gap) stays a scalar under its
    // literal name; guessing would invent components.
    cgns_mll_read_solutions(fn, B, Z, npoints, mesh, rPath);

    return mesh;
}

#endif  // MESHIOPLUSPLUS_HAS_CGNSLIB

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
