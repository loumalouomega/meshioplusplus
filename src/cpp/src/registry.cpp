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

/**
 * @file registry.cpp
 * @brief The shared format-dispatch tables (see registry.hpp). Bodies hoisted
 *        verbatim from `bindings/wasm/js_bindings.cpp`, extended with the
 *        HDF5/netCDF-conditional entries native (non-WASM) builds can serve.
 */

// System includes
#include <unordered_map>

// Project includes
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/formats/triangle.hpp"
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vti.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/xdmf.hpp"

namespace meshioplusplus {

const std::map<std::string, ReadFn>& registry_readers() {
    static const std::map<std::string, ReadFn> m = {
        {"abaqus", meshioplusplus::read_abaqus},
        {"ansys", meshioplusplus::read_ansys},
        {"avsucd", meshioplusplus::read_avsucd},
        {"dolfin", meshioplusplus::read_dolfin},
        {"ensight", meshioplusplus::read_ensight},
        {"flac3d", meshioplusplus::read_flac3d},
        {"dex", meshioplusplus::read_dex},
        {"flux", meshioplusplus::read_flux},
        {"freefem", meshioplusplus::read_freefem},
        // UNGUARDED, unlike gid's writer entry: gidpost is write-only, so
        // reading needs none of it. read_gid's real dependencies are per
        // flavour (ascii: none, binary: zlib, hdf5: HDF5), which makes `gid`
        // readable in strictly more build configurations than it is writable
        // -- the release CLI binaries and Windows wheels build zlib-off and
        // cannot write GiD at all, yet read the ascii flavour fine. For the
        // same reason gid must NOT appear in registry_compiled_out(): there is
        // no missing dependency to name.
        {"gid", [](const std::string& p) { return meshioplusplus::read_gid(p); }},
        {"gmsh", [](const std::string& path) { return meshioplusplus::read_gmsh(path); }},
        {"ip", meshioplusplus::read_ip},
        // mdpa now has read overloads (ReadOptions / MdpaInfo), so the plain
        // function pointer no longer converts to ReadFn -- the hazard noted
        // for unv/med below. The MdpaInfo is dropped here, like MedInfo.
        {"mdpa", [](const std::string& path) { return meshioplusplus::read_mdpa(path); }},
        {"medit", meshioplusplus::read_medit_ascii},
        {"mff", meshioplusplus::read_mff},
        {"mfm", meshioplusplus::read_mfm},
        {"mphtxt", meshioplusplus::read_mphtxt},
        {"nastran", meshioplusplus::read_nastran},
        {"netgen", meshioplusplus::read_netgen},
        {"obj", meshioplusplus::read_obj},
        {"off", meshioplusplus::read_off},
        {"permas", meshioplusplus::read_permas},
        {"ply", meshioplusplus::read_ply},
        {"stl", meshioplusplus::read_stl},
        {"su2", meshioplusplus::read_su2},
        {"tecplot", meshioplusplus::read_tecplot},
        {"tetgen", meshioplusplus::read_tetgen},
        {"triangle", meshioplusplus::read_triangle},
        {"ugrid", meshioplusplus::read_ugrid},
        {"unv", [](const std::string& path) { return meshioplusplus::read_unv(path); }},
        {"vti", [](const std::string& path) { return meshioplusplus::read_vti(path); }},
        {"vtk", meshioplusplus::read_vtk},
        // vti/vtp/vtu take a trailing defaulted ReadOptions, so the function
        // pointers no longer convert to ReadFn -- wrapped like unv/med below.
        {"vtp", [](const std::string& path) { return meshioplusplus::read_vtp(path); }},
        {"vtu", [](const std::string& path) { return meshioplusplus::read_vtu(path); }},
        {"wkt", meshioplusplus::read_wkt},
        {"xdmf", [](const std::string& path) { return meshioplusplus::read_xdmf(path); }},
        // Side-channel info (point_sets/cell_sets, cell-tag family names) is
        // not carried by the flat bindings -- v1 limitation, see doc/wasm.md
        // and doc/c_api.md.
        {"ansysinp",
         [](const std::string& path) {
             meshioplusplus::AnsysInfo info;
             return meshioplusplus::read_ansysinp(path, info);
         }},
        {"openfoam",
         [](const std::string& path) {
             meshioplusplus::OpenFoamInfo info;
             return meshioplusplus::read_openfoam(path, info);
         }},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        {"cgns", meshioplusplus::read_cgns},
        {"h5m", meshioplusplus::read_h5m},
        {"hmf", meshioplusplus::read_hmf},
        {"med",
         [](const std::string& path) {
             // The family-id maps/link names/mesh metadata in MedInfo are
             // still dropped here, but group *names* are not lost: read_med
             // attaches them as named regions directly on the Mesh (see
             // med_attach_point_regions/med_attach_cell_regions in med.cpp),
             // so they reach WASM/C API/Fortran through this path too.
             meshioplusplus::MedInfo info;
             return meshioplusplus::read_med(path, info);
         }},
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        // Explicit lambda, not `&read_exodus`: the overload set now also holds
        // the ExodusInfo form, and taking its address would be ambiguous. The
        // provenance strings are dropped here, as MedInfo's are below.
        {"exodus", [](const std::string& path) { return meshioplusplus::read_exodus(path); }},
#endif
    };
    return m;
}

const std::map<std::string, WriteFn>& registry_writers() {
    static const std::map<std::string, WriteFn> m = {
        {"abaqus", meshioplusplus::write_abaqus},
        {"ansys", [](const std::string& p,
                     const Mesh& mm) { meshioplusplus::write_ansys(p, mm, /*binary=*/true); }},
        {"avsucd", meshioplusplus::write_avsucd},
        {"dolfin", meshioplusplus::write_dolfin},
        {"ensight", [](const std::string& p,
                       const Mesh& mm) { meshioplusplus::write_ensight(p, mm, /*binary=*/true); }},
        {"flac3d",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_flac3d(p, mm, ".16e", /*binary=*/false);
         }},
        {"dex", meshioplusplus::write_dex},
        {"flux", meshioplusplus::write_flux},
        {"freefem", meshioplusplus::write_freefem},
        // Write-only (gidpost has no read functions at all -- the reader is a
        // documented follow-up, doc/roadmap.md section 1). GidMode::Auto
        // infers the flavour (ascii/binary/hdf5) from the path's extension.
        {"gid", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_gid(p, mm); }},
        {"gmsh", [](const std::string& p,
                    const Mesh& mm) { meshioplusplus::write_gmsh41(p, mm, /*binary=*/true); }},
        // Distinct from "gmsh" (4.1): 2.2 stores each element's physical tag
        // directly, so it is the version that round-trips named Cell region
        // MEMBERSHIP, not just names -- write_gmsh22 already synthesizes
        // gmsh:physical from Cell regions when the mesh has none of its own.
        // Was reachable only from Python (gmsh22_write) until this entry; the
        // flat bindings (WASM/C API/Fortran) had no way to select it at all.
        {"gmsh22", [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_gmsh22(p, mm, /*binary=*/true); }},
        {"ip", meshioplusplus::write_ip},
        {"mdpa", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_mdpa(p, mm); }},
        {"medit", meshioplusplus::write_medit_ascii},
        {"mff", meshioplusplus::write_mff},
        {"mfm",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_mfm(p, mm, ".16e"); }},
        {"mphtxt", meshioplusplus::write_mphtxt},
        {"nastran", meshioplusplus::write_nastran},
        {"netgen",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_netgen(p, mm, ".16e"); }},
        {"obj", meshioplusplus::write_obj},
        {"off", meshioplusplus::write_off},
        {"permas", meshioplusplus::write_permas},
        // ply/stl default to skin=true (matching the Python shims): a volume
        // mesh writes its extracted boundary skin instead of dropping the
        // volume cells.
        {"ply",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_ply(p, mm, /*binary=*/true, /*skin=*/true);
         }},
        {"stl",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_stl(p, mm, /*binary=*/false, /*skin=*/true);
         }},
        {"su2", meshioplusplus::write_su2},
        // svg/tikz are write-only 2D-visualization formats; the flat bindings
        // emit them with the fixed default styling (per-call overrides are out
        // of scope for v1, per registry.hpp).
        {"svg", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_svg(p, mm); }},
        {"tikz", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_tikz(p, mm); }},
        {"tecplot", meshioplusplus::write_tecplot},
        {"tetgen", meshioplusplus::write_tetgen},
        {"triangle", meshioplusplus::write_triangle},
        {"ugrid", meshioplusplus::write_ugrid},
        {"unv", [](const std::string& p, const Mesh& mm) { meshioplusplus::write_unv(p, mm); }},
        {"vti",
         [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
             meshioplusplus::write_vti(p, mm, /*binary=*/true, /*zlib=*/true);
#else
             meshioplusplus::write_vti(p, mm, /*binary=*/true, /*zlib=*/false);
#endif
         }},
        {"vtk",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_vtk(p, mm, /*binary=*/true, /*v51=*/true);
         }},
        // The zlib codec follows the build, exactly like the xdmf entry below:
        // a hardcoded zlib=true is a WriteError on a -DMESHIOPLUSPLUS_WITH_ZLIB=OFF
        // build (comment above already documented "falls back to Python" as
        // the assumption -- true for every binding with a Python shim to fall
        // back to, but not for a caller reaching this table directly, e.g.
        // WASM/C API/Fortran or the sequence engine's per-step writes).
        {"vtp",
         [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
             meshioplusplus::write_vtp(p, mm, /*binary=*/true, /*zlib=*/true);
#else
             meshioplusplus::write_vtp(p, mm, /*binary=*/true, /*zlib=*/false);
#endif
         }},
        {"vtu",
         [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
             meshioplusplus::write_vtu(p, mm, /*binary=*/true, /*zlib=*/true);
#else
             meshioplusplus::write_vtu(p, mm, /*binary=*/true, /*zlib=*/false);
#endif
         }},
        {"wkt", meshioplusplus::write_wkt},
        // XDMF's heavy-data format follows the build: HDF companion file when
        // HDF5 is available (the Python writer's default), inline XML text
        // otherwise (the only always-available option; what WASM ships).
        {"xdmf",
         [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
             meshioplusplus::write_xdmf(p, mm, "HDF");
#else
             meshioplusplus::write_xdmf(p, mm, "XML");
#endif
         }},
        {"ansysinp",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::AnsysInfo info;  // no point_sets/cell_sets side channel in v1
             meshioplusplus::write_ansysinp(p, mm, info);
         }},
        {"openfoam",
         [](const std::string& p, const Mesh& mm) {
             // No patch-name/type side channel here, so this yields a single
             // `defaultFaces` patch -- a valid, loadable case (see
             // write_openfoam). The flat bindings get that rather than an error.
             meshioplusplus::OpenFoamInfo info;
             meshioplusplus::write_openfoam(p, mm, info);
         }},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        {"cgns", [](const std::string& p,
                    const Mesh& mm) { meshioplusplus::write_cgns(p, mm, /*gzip_level=*/4); }},
        {"h5m",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_h5m(p, mm, /*add_global_ids=*/true, /*gzip_level=*/4);
         }},
        {"hmf", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_hmf(p, mm, /*gzip_level=*/4); }},
        {"med",
         [](const std::string& p, const Mesh& mm) {
             // No point_tags/cell_tags to hand over here (they live in the
             // dropped MedInfo), but write_med synthesizes them from any
             // Point/Cell regions the mesh carries (see
             // med_point_regions_to_tags/med_cell_regions_to_tags in
             // med.cpp), so a mesh converted from e.g. Abaqus through this
             // registry path still carries its named groups into the file.
             meshioplusplus::MedInfo info;
             meshioplusplus::write_med(p, mm, info);
         }},
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        {"exodus", meshioplusplus::write_exodus},
#endif
    };
    return m;
}

// Extension -> canonical format key for the non-ambiguous cases; `.msh`
// defaults to gmsh and `.inp` to abaqus (matching this repo's own import
// order in src/python/meshioplusplus/__init__.py). Pass an explicit `format` to
// select ansys/freefem (.msh) or ansysinp (.inp) instead. Optional-dependency
// extensions are mapped even in builds where the format is compiled out, so
// the resulting error names the missing dependency (registry_compiled_out())
// rather than claiming the extension is unknown.
const std::map<std::string, std::string>& registry_extension_defaults() {
    static const std::map<std::string, std::string> m = {
        {".inp", "abaqus"},
        {".avs", "avsucd"},
        {".xml", "dolfin"},
        {".f3grid", "flac3d"},
        {".case", "ensight"},
        {".geo", "ensight"},
        {".dex", "dex"},
        {".ip", "ip"},
        {".mff", "mff"},
        {".pf3", "flux"},
        {".mdpa", "mdpa"},
        {".mesh", "medit"},
        {".mfm", "mfm"},
        {".mphtxt", "mphtxt"},
        {".bdf", "nastran"},
        {".nas", "nastran"},
        {".fem", "nastran"},
        {".vol", "netgen"},
        {".obj", "obj"},
        // OpenFOAM: the `.foam` marker file. A case *directory* has no
        // extension at all, so that form still needs an explicit format.
        {".foam", "openfoam"},
        {".off", "off"},
        {".post", "permas"},
        {".dato", "permas"},
        // GiD postprocess. All four are compound extensions, and
        // resolve_format() tries the longest suffix first, so ".post.msh"
        // resolves to "gid" rather than falling through to ".msh" -> gmsh.
        // Registered even when the build has no gidpost: write_gid() then
        // throws naming the missing build flags, which is strictly better
        // than ".post.msh" silently resolving to another format.
        {".post.msh", "gid"},
        {".post.res", "gid"},
        {".post.bin", "gid"},
        {".post.h5", "gid"},
        {".ply", "ply"},
        {".stl", "stl"},
        {".su2", "su2"},
        {".svg", "svg"},
        {".tikz", "tikz"},
        // .node/.ele stay with tetgen for backward compatibility (3D pairs
        // are the common case for the flat bindings); reading Triangle 2D
        // .node/.ele files there needs an explicit format="triangle". Only
        // .poly defaults to triangle.
        {".dat", "tecplot"},
        {".tec", "tecplot"},
        {".ele", "tetgen"},
        {".node", "tetgen"},
        {".poly", "triangle"},
        {".ugrid", "ugrid"},
        {".unv", "unv"},
        {".vti", "vti"},
        {".vtk", "vtk"},
        {".vtp", "vtp"},
        {".vtu", "vtu"},
        {".wkt", "wkt"},
        {".xdmf", "xdmf"},
        {".xmf", "xdmf"},
        {".msh", "gmsh"},
        {".cgns", "cgns"},
        {".h5m", "h5m"},
        {".hmf", "hmf"},
        {".med", "med"},
        {".e", "exodus"},
        {".exo", "exodus"},
        {".ex2", "exodus"},
    };
    return m;
}

namespace {

// Longest compound extension first: ".post.msh" must win over ".msh" (gmsh),
// or the GiD writer is unreachable by path alone. Walking from the FIRST dot
// of the basename forward yields candidates in strictly decreasing length,
// so the first hit is the longest match. The basename strip matters: a
// directory component with a dot ("/home/.config/m.vtu") would otherwise
// produce a nonsense first candidate. Behaviour-preserving for every
// extension registered before the compound ones existed -- every one of them
// is single-dot, so at most one candidate can ever match for those paths.
std::string basename_of(const std::string& rPath) {
    auto pos = rPath.find_last_of("/\\");
    return pos == std::string::npos ? rPath : rPath.substr(pos + 1);
}

}  // namespace

std::string resolve_format(const std::string& rPath, const std::string& rFormat) {
    if (!rFormat.empty())
        return rFormat;
    const auto& defaults = registry_extension_defaults();
    const std::string base = basename_of(rPath);
    for (std::size_t pos = base.find('.'); pos != std::string::npos;
         pos = base.find('.', pos + 1)) {
        auto it = defaults.find(base.substr(pos));
        if (it != defaults.end())
            return it->second;
    }
    throw meshioplusplus::ReadError("meshio++: cannot infer format from '" + rPath +
                                    "' -- pass an explicit format argument");
}

const std::unordered_map<std::string, ReadExFn>& registry_readers_ex() {
    // Sparse by design -- populated per format as native selective-read support
    // lands. An absent format falls back to a full read in registry_read().
    static const std::unordered_map<std::string, ReadExFn> m = {
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        // Exodus honours mTimeStep rather than the narrowing options -- being
        // "options-aware" is not one capability but several, and a time series
        // is the one this format has. IWYU pragma: keep
        {"exodus", [](const std::string& path,
                      const ReadOptions& opts) { return meshioplusplus::read_exodus(path, opts); }},
#endif
        // A lambda, not `&read_gmsh`: the GmshInfo overload makes the bare name
        // ambiguous (the exodus/mdpa story again). The info is dropped here, so
        // the flat bindings see no `$Entities` bounding entities.
        {"gmsh", [](const std::string& path,
                    const ReadOptions& opts) { return meshioplusplus::read_gmsh(path, opts); }},
        // mdpa honours mLenient rather than the narrowing options -- the same
        // "options-aware is several capabilities, not one" note as exodus. This
        // entry is what makes `--lenient` reach the C API, Fortran, Julia, R,
        // WASM and the native CLI with no per-binding code.
        {"mdpa", [](const std::string& path,
                    const ReadOptions& opts) { return meshioplusplus::read_mdpa(path, opts); }},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        // MED honours `mLenient` (skip/report the enhanced `CHA` constructs
        // instead of deferring the whole file to Python) and `mTimeStep`
        // (select one step of a multi-step field) -- neither is a narrowing
        // option, the same "options-aware is several capabilities, not one"
        // note as exodus and mdpa. The `MedInfo` is dropped here exactly as
        // the plain reader entry drops it, so the flat bindings get the mesh
        // and its representable fields but not `mFieldUnits`/`mStepMeta`; that
        // is a documented gap, not a silent loss, and it is strictly more than
        // the hard failure they used to get. IWYU pragma: keep
        {"med",
         [](const std::string& path, const ReadOptions& opts) {
             meshioplusplus::MedInfo info;
             return meshioplusplus::read_med(path, info, opts);
         }},
#endif
        {"gid", meshioplusplus::read_gid},
        {"vti", meshioplusplus::read_vti},
        {"vtp", meshioplusplus::read_vtp},
        {"vtu", meshioplusplus::read_vtu},
        {"xdmf", meshioplusplus::read_xdmf},
    };
    return m;
}

const std::unordered_map<std::string, MetadataFn>& registry_metadata_readers() {
    static const std::unordered_map<std::string, MetadataFn> m = {
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        {"exodus", meshioplusplus::read_exodus_metadata},
#endif
        {"gmsh", meshioplusplus::read_gmsh_metadata},
        {"gid", meshioplusplus::read_gid_metadata},
        {"vti", meshioplusplus::read_vti_metadata},
        {"vtp", meshioplusplus::read_vtp_metadata},
        {"vtu", meshioplusplus::read_vtu_metadata},
        {"xdmf", meshioplusplus::read_xdmf_metadata},
    };
    return m;
}

bool registry_reader_supports_options(const std::string& rFormat) {
    return registry_readers_ex().count(rFormat) > 0;
}

namespace {

/** @brief The reader for @p rFormat, or a ReadError naming why it is missing. */
const ReadFn& registry_full_reader(const std::string& rFormat) {
    auto it = registry_readers().find(rFormat);
    if (it == registry_readers().end()) {
        const char* dep = registry_compiled_out(rFormat);
        throw meshioplusplus::ReadError(
            "meshio++: unknown or unsupported format '" + rFormat + "'" +
            (dep ? std::string(" (this build has no ") + dep + " support)" : std::string()));
    }
    return it->second;
}

}  // namespace

Mesh registry_read(const std::string& rPath, const std::string& rFormat,
                   const ReadOptions& rOptions) {
    auto it = registry_readers_ex().find(rFormat);
    if (it != registry_readers_ex().end())
        return it->second(rPath, rOptions);
    // No native selective path: a full read is still the correct answer.
    return registry_full_reader(rFormat)(rPath);
}

MeshMetadata registry_read_metadata(const std::string& rPath, const std::string& rFormat,
                                    const ReadOptions& rOptions) {
    // Best-effort, on every path: the block lives in the file's bytes, so
    // unlike everything else in a summary it cannot come from
    // `metadata_from_mesh`. A reader that fills it natively (exodus, whose
    // block is a netCDF attribute rather than head bytes) wins -- hence the
    // "only if still empty" test below rather than an unconditional overwrite.
    auto fill_provenance = [&rPath](MeshMetadata& rMeta) {
        if (!rMeta.mProvenance.empty())
            return;
        detail::ProvenanceReadResult found = detail::read_provenance_lines(rPath);
        rMeta.mProvenance = std::move(found.mLines);
        rMeta.mProvenanceRecognised = found.mRecognised;
    };

    auto it = registry_metadata_readers().find(rFormat);
    if (it != registry_metadata_readers().end()) {
        try {
            MeshMetadata meta = it->second(rPath, rOptions);
            meta.mFormat = rFormat;
            fill_provenance(meta);
            return meta;
        } catch (const meshioplusplus::ReadError&) {
            // A native summary can legitimately decline a construct it cannot
            // describe cheaply but the full reader handles fine (XDMF `Mixed`
            // topology is the motivating case: per-block counts are only
            // knowable after reading the topology array). Declining must cost a
            // slower answer, not a failed one -- so fall through to the full
            // read. A genuinely unreadable file still throws below.
        }
    }
    // The one fallback for every format lacking a native metadata path -- read
    // it whole, summarize, and say so rather than implying it was cheap.
    MeshMetadata meta = metadata_from_mesh(registry_full_reader(rFormat)(rPath));
    meta.mFellBackToFullRead = true;
    meta.mFormat = rFormat;
    fill_provenance(meta);
    return meta;
}

const char* registry_compiled_out(const std::string& rFormat) {
#ifndef MESHIOPLUSPLUS_HAS_HDF5
    if (rFormat == "cgns" || rFormat == "h5m" || rFormat == "hmf" || rFormat == "med")
        return "HDF5";
#endif
#ifndef MESHIOPLUSPLUS_HAS_NETCDF
    if (rFormat == "exodus")
        return "netCDF";
#endif
    // Deliberately NOT an arm for cgns/cgnslib: the format is fully readable
    // and writable without it (the hand-rolled ADF-over-HDF5 path), so
    // reporting it "compiled out" would be wrong. The cgnslib-only
    // capabilities -- ADF containers and NGON_n/NFACE_n -- report themselves
    // by name from read_cgns_mll instead.
    return nullptr;
}

}  // namespace meshioplusplus
