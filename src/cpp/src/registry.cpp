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
#include <cctype>
#include <filesystem>
#include <ios>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/compressed_input.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/read_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/abaqus_fil.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/lsdyna.hpp"
#include "meshioplusplus/formats/lsdyna_binout.hpp"
#include "meshioplusplus/formats/lsdyna_d3plot.hpp"
#include "meshioplusplus/formats/code_aster.hpp"
#include "meshioplusplus/formats/patran.hpp"
#include "meshioplusplus/formats/elmer.hpp"
#include "meshioplusplus/formats/febio.hpp"
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/formats/libmesh.hpp"
#include "meshioplusplus/formats/z88.hpp"
#include "meshioplusplus/formats/radioss.hpp"
#include "meshioplusplus/formats/radioss_anim.hpp"
#include "meshioplusplus/formats/radioss_th.hpp"
#include "meshioplusplus/formats/xplt.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ansys_rst.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/adcirc.hpp"
#include "meshioplusplus/formats/gifti.hpp"
#include "meshioplusplus/formats/freesurfer.hpp"
#include "meshioplusplus/formats/selafin.hpp"
#include "meshioplusplus/formats/dxf.hpp"
#include "meshioplusplus/formats/plot3d.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/frd.hpp"
#include "meshioplusplus/formats/ugrid_nc.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/marc.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/nastran_op2.hpp"
#include "meshioplusplus/formats/nastran_h5.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/pcd.hpp"
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
#include "meshioplusplus/formats/vts.hpp"
#include "meshioplusplus/formats/vtr.hpp"
#include "meshioplusplus/formats/pvd.hpp"
#include "meshioplusplus/formats/pvtp.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/formats/vtx.hpp"
#include "meshioplusplus/formats/szplt.hpp"
#include "meshioplusplus/formats/vtm.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/formats/neuroglancer.hpp"
#include "meshioplusplus/formats/xyz.hpp"

namespace meshioplusplus {

namespace {

/**
 * @brief Wraps every reader so a parser's std:: exception leaves as ReadError.
 *
 * Every consumer of the registry -- the native CLI, the C API, WASM and the
 * fuzz targets -- then sees one exception for "not a file this reader reads"
 * (detail/read_guard.hpp explains why that matters to the Python fallback).
 */
// The codec a path's compression suffix names, unless a format registers the
// whole compound extension itself (Netgen's `.vol.gz`, libMesh's `.xda.gz`):
// that one decompresses its own files. See detail/compressed_input.hpp.
detail::Compression registry_compression(const std::string& rPath) {
    const detail::Compression kind = detail::compression_from_suffix(rPath);
    if (kind == detail::Compression::None)
        return kind;
    std::string base = std::filesystem::path(rPath).filename().string();
    for (char& c : base)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    const auto& defaults = registry_extension_defaults();
    const std::string own = detail::compression_suffix(kind);
    for (std::size_t pos = base.find('.'); pos != std::string::npos;
         pos = base.find('.', pos + 1)) {
        const std::string suffix = base.substr(pos);
        if (suffix != own && defaults.count(suffix))
            return detail::Compression::None;
    }
    return kind;
}

// Formats whose file names other files beside it (or is a directory): a
// decompressed copy in a temporary directory would lose them, so a compressed
// one is refused by name rather than read half.
bool registry_is_multi_file(const std::string& rFormat) {
    static const char* const kMulti[] = {"tetgen", "triangle", "ensight", "gid", "elmer",
                                         "pvd",    "pvtu",     "pvtp",    "vtx", "openfoam"};
    for (const char* name : kMulti)
        if (rFormat == name)
            return true;
    return false;
}

void registry_refuse_multi_file(const std::string& rFormat, const std::string& rPath) {
    if (registry_is_multi_file(rFormat))
        throw meshioplusplus::ReadError(
            "meshio++: '" + rPath + "' is a compressed " + rFormat +
            " file, whose data spans several files; decompress it first");
}

std::map<std::string, ReadFn> registry_guard_readers(std::map<std::string, ReadFn> raw) {
    for (auto& [name, fn] : raw) {
        fn = [format = name, inner = std::move(fn)](const std::string& rPath) {
            // A compressed file is read through a decompressed copy, so every
            // flat binding and both CLIs read `cube.msh.gz` with no code of
            // their own (roadmap 7.1.1).
            const detail::Compression kind = registry_compression(rPath);
            if (kind != detail::Compression::None) {
                registry_refuse_multi_file(format, rPath);
                const detail::DecompressedFile plain(rPath, kind);
                return detail::guarded_read(format.c_str(), [&] { return inner(plain.Path()); });
            }
            return detail::guarded_read(format.c_str(), [&] { return inner(rPath); });
        };
    }
    return raw;
}

}  // namespace

namespace {

// Whether the VTK XML writers compress their arrays. The zlib codec follows
// the build, exactly like the xdmf entry's heavy-data format: a hardcoded
// zlib=true is a WriteError on a -DMESHIOPLUSPLUS_WITH_ZLIB=OFF build, and
// "falls back to Python" holds for every binding with a Python shim, but not
// for a caller reaching this table directly, e.g. WASM/C API/Fortran or the
// sequence engine's per-step writes.
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
constexpr bool kRegistryZlib = true;
#else
constexpr bool kRegistryZlib = false;
#endif

/**
 * @brief One format's entries in the registry: what reads it, what writes it,
 * its options-aware reader and its metadata reader.
 *
 * An empty function means the format does not offer that entry (gmsh22 only
 * writes, iges only reads, most formats have no options-aware reader).
 * registry_readers(), registry_writers(), registry_readers_ex() and
 * registry_metadata_readers() are built from these rows; the extension table,
 * many extensions to one format, stays its own.
 */
struct RegistryRow {
    const char* mName;
    ReadFn mRead;
    WriteFn mWrite;
    ReadExFn mReadEx;
    MetadataFn mMetadata;
};

const std::vector<RegistryRow>& registry_rows() {
    static const std::vector<RegistryRow> rows = {
        {.mName = "abaqus",
         .mRead = meshioplusplus::read_abaqus,
         .mWrite = meshioplusplus::write_abaqus},
        {.mName = "abaqus_fil",
         .mRead = [](const std::string& path) { return meshioplusplus::read_abaqus_fil(path); },
         .mWrite = [](const std::string& p,
                      const Mesh& m) { meshioplusplus::write_abaqus_fil(p, m); },
         // Abaqus .fil honours mTimeStep (its steps are the increments) and the
         // narrowing options.
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_abaqus_fil(path, opts);
             },
         .mMetadata = meshioplusplus::read_abaqus_fil_metadata},
        {.mName = "lsdyna",
         .mRead = meshioplusplus::read_lsdyna,
         .mWrite = meshioplusplus::write_lsdyna},
        // The d3plot family's base file is named `d3plot`: resolve_format
        // matches the basename, sniff_format the control block.
        {.mName = "lsdyna_d3plot",
         .mRead = [](const std::string& path) { return meshioplusplus::read_lsdyna_d3plot(path); },
         .mWrite = [](const std::string& rPath,
                      const Mesh& rMesh) { meshioplusplus::write_lsdyna_d3plot(rPath, rMesh); },
         // LS-DYNA d3plot honours mTimeStep (its steps are the family's states)
         // and the narrowing options.
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_lsdyna_d3plot(path, opts);
             },
         .mMetadata = meshioplusplus::read_lsdyna_d3plot_metadata},
        // LS-DYNA's binary output database: `binout`, found by name or header.
        {.mName = "lsdyna_binout",
         .mRead = [](const std::string& path) { return meshioplusplus::read_lsdyna_binout(path); },
         // No extension (`binout`): by name or explicit format only.
         .mWrite = [](const std::string& p,
                      const Mesh& m) { meshioplusplus::write_lsdyna_binout(p, m); },
         // LS-DYNA binout honours mTimeStep (its steps are nodout's outputs)
         // and the narrowing options.
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_lsdyna_binout(path, opts);
             },
         .mMetadata = meshioplusplus::read_lsdyna_binout_metadata},
        {.mName = "code_aster",
         .mRead = meshioplusplus::read_code_aster,
         .mWrite = meshioplusplus::write_code_aster},
        // A lambda: read_patran is overloaded (result files).
        {.mName = "patran",
         .mRead = [](const std::string& path) { return meshioplusplus::read_patran(path); },
         .mWrite = meshioplusplus::write_patran},
        {.mName = "femap",
         .mRead = [](const std::string& path) { return meshioplusplus::read_femap(path); },
         .mWrite = meshioplusplus::write_femap,
         // Femap honours mTimeStep (its steps are the 450 output sets) and the
         // data narrowing options.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_femap(path, opts); },
         .mMetadata = meshioplusplus::read_femap_metadata},
        {.mName = "libmesh",
         .mRead = meshioplusplus::read_libmesh,
         .mWrite = meshioplusplus::write_libmesh},
        {.mName = "radioss",
         .mRead = meshioplusplus::read_radioss,
         .mWrite = [](const std::string& p,
                      const Mesh& m) { meshioplusplus::write_radioss(p, m); }},
        // Found by name (`<run>A001`...) or magic, never an extension.
        {.mName = "radioss_anim",
         .mRead = meshioplusplus::read_radioss_anim,
         .mWrite = meshioplusplus::write_radioss_anim,
         .mMetadata = meshioplusplus::read_radioss_anim_metadata},
        // OpenRadioss time-history files: `<run>T01`, found by name or header.
        {.mName = "radioss_th",
         .mRead = [](const std::string& path) { return meshioplusplus::read_radioss_th(path); },
         // No extension (`<run>T01`): by name or explicit format only.
         .mWrite = meshioplusplus::write_radioss_th,
         // Radioss time history honours mTimeStep (its steps are the outputs)
         // and the narrowing options.
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_radioss_th(path, opts);
             },
         .mMetadata = meshioplusplus::read_radioss_th_metadata},
        // A .dat file is Tecplot's unless it opens as a Marc deck (resolve_format).
        {.mName = "marc",
         .mRead = meshioplusplus::read_marc,
         // By name only: ".dat" writes Tecplot (resolve_write_format).
         .mWrite = meshioplusplus::write_marc},
        {.mName = "marc_t19",
         .mRead = [](const std::string& path) { return meshioplusplus::read_marc_t19(path); },
         .mWrite = meshioplusplus::write_marc_t19,
         // Marc .t19: mTimeStep picks the increment.
         .mReadEx = meshioplusplus::read_marc_t19,
         .mMetadata = meshioplusplus::read_marc_t19_metadata},
        // Fixed file names (z88i1.txt ...): resolve_format matches the basename.
        {.mName = "z88",
         .mRead = [](const std::string& path) { return meshioplusplus::read_z88(path); },
         .mWrite = [](const std::string& p, const Mesh& m) { meshioplusplus::write_z88(p, m); }},
        // A directory, not a file: no extension maps to it; sniff_format finds it.
        {.mName = "elmer",
         .mRead = [](const std::string& path) { return meshioplusplus::read_elmer(path); },
         .mWrite = [](const std::string& p, const Mesh& m) { meshioplusplus::write_elmer(p, m); },
         // Elmer honours mPiece/mPieceSet (one part of a partitioned mesh) and
         // mLenient (skip element types with no meshio++ cell type).
         .mReadEx = meshioplusplus::read_elmer},
        {.mName = "febio",
         .mRead = [](const std::string& path) { return meshioplusplus::read_febio(path); },
         .mWrite = meshioplusplus::write_febio,
         // FEBio .feb honours mLenient (downgrade tet5/tet15 to tetra/tetra10).
         .mReadEx = meshioplusplus::read_febio},
        {.mName = "xplt",
         .mRead = [](const std::string& path) { return meshioplusplus::read_xplt(path); },
         .mWrite = [](const std::string& rPath,
                      const Mesh& rMesh) { meshioplusplus::write_xplt(rPath, rMesh); },
         // FEBio .xplt honours mTimeStep (one state), the narrowing options and
         // mLenient (downgrade tet5/tet15 domains).
         .mReadEx = meshioplusplus::read_xplt,
         .mMetadata = meshioplusplus::read_xplt_metadata},
        {.mName = "ansys_rst",
         .mRead = [](const std::string& path) { return meshioplusplus::read_ansys_rst(path); },
         .mWrite = meshioplusplus::write_ansys_rst,
         // Ansys .rst/.rth: mTimeStep picks the result set, like .xplt.
         .mReadEx = meshioplusplus::read_ansys_rst,
         .mMetadata = meshioplusplus::read_ansys_rst_metadata},
        {.mName = "ansys_rst_cyclic",
         .mRead =
             [](const std::string& path) { return meshioplusplus::read_ansys_rst_cyclic(path); },
         // Its full-rotor reading (a static cyclic model): no extension, by name only.
         .mReadEx = meshioplusplus::read_ansys_rst_cyclic,
         .mMetadata = meshioplusplus::read_ansys_rst_cyclic_metadata},
        // Read-only, and a lambda for the same reason as ensight's: overloaded.
        {.mName = "frd",
         .mRead = [](const std::string& path) { return meshioplusplus::read_frd(path); },
         .mWrite = [](const std::string& p, const Mesh& m) { meshioplusplus::write_frd(p, m); },
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_frd(path, opts); },
         .mMetadata = meshioplusplus::read_frd_metadata},
        {.mName = "ansys",
         .mRead = meshioplusplus::read_ansys,
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_ansys(p, mm, /*binary=*/true); }},
        {.mName = "avsucd",
         .mRead = meshioplusplus::read_avsucd,
         .mWrite = meshioplusplus::write_avsucd},
        {.mName = "dolfin",
         .mRead = meshioplusplus::read_dolfin,
         .mWrite = meshioplusplus::write_dolfin},
        // A lambda, not `&read_ensight`: the ReadOptions overload makes the
        // bare name ambiguous (the exodus/mdpa/med/cgns/tecplot story again).
        {.mName = "ensight",
         .mRead = [](const std::string& path) { return meshioplusplus::read_ensight(path); },
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_ensight(p, mm, /*binary=*/true); },
         // EnSight honours mTimeStep AND the narrowing options -- a .case
         // file's VARIABLE entries are only ever read here, never by the
         // plain overload. IWYU pragma: keep
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_ensight(path, opts);
             },
         .mMetadata = meshioplusplus::read_ensight_metadata},
        {.mName = "flac3d",
         .mRead = meshioplusplus::read_flac3d,
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_flac3d(p, mm, ".16e", /*binary=*/false);
             }},
        {.mName = "dex", .mRead = meshioplusplus::read_dex, .mWrite = meshioplusplus::write_dex},
        {.mName = "flux", .mRead = meshioplusplus::read_flux, .mWrite = meshioplusplus::write_flux},
        {.mName = "freefem",
         .mRead = meshioplusplus::read_freefem,
         .mWrite = meshioplusplus::write_freefem},
        // Kratos CAD JSON: a lambda, read_cad_json is overloaded (BrepInfo).
        {.mName = "cad_json",
         .mRead = [](const std::string& p) { return meshioplusplus::read_cad_json(p); },
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_cad_json(p, mm); }},
        // The reader is UNGUARDED, unlike the writer: gidpost is write-only, so
        // reading needs none of it. read_gid's real dependencies are per
        // flavour (ascii: none, binary: zlib, hdf5: HDF5), which makes `gid`
        // readable in strictly more build configurations than it is writable
        // -- the release CLI binaries and Windows wheels build zlib-off and
        // cannot write GiD at all, yet read the ascii flavour fine. For the
        // same reason gid must NOT appear in registry_compiled_out(): there is
        // no missing dependency to name.
        {.mName = "gid",
         .mRead = [](const std::string& p) { return meshioplusplus::read_gid(p); },
         // gidpost itself has no read functions at all -- meshio++'s own
         // read_gid (the mRead entry above) is a hand-rolled reader that
         // does not depend on it, which is why gid is readable in strictly
         // more build configurations than it is writable. GidMode::Auto
         // infers the flavour (ascii/binary/hdf5) from the path's extension.
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_gid(p, mm); },
         .mReadEx = meshioplusplus::read_gid,
         .mMetadata = meshioplusplus::read_gid_metadata},
        // IGES and STEP read through the optional OpenCASCADE backend
        // (roadmap 6.2.1); registered either way, so a build without it
        // throws Unsupported naming the flag instead of "unknown format".
        {.mName = "iges",
         .mRead = [](const std::string& p) { return meshioplusplus::read_iges(p); }},
        {.mName = "gmsh",
         .mRead = [](const std::string& path) { return meshioplusplus::read_gmsh(path); },
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_gmsh41(p, mm, /*binary=*/true); },
         // A lambda, not `&read_gmsh`: the GmshInfo overload makes the bare name
         // ambiguous (the exodus/mdpa story again). The info is dropped here, so
         // the flat bindings see no `$Entities` bounding entities.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_gmsh(path, opts); },
         .mMetadata = meshioplusplus::read_gmsh_metadata},
        {.mName = "ip", .mRead = meshioplusplus::read_ip, .mWrite = meshioplusplus::write_ip},
        // mdpa now has read overloads (ReadOptions / MdpaInfo), so the plain
        // function pointer no longer converts to ReadFn -- the hazard noted
        // for unv/med. The MdpaInfo is dropped here, like MedInfo.
        {.mName = "mdpa",
         .mRead = [](const std::string& path) { return meshioplusplus::read_mdpa(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_mdpa(p, mm); },
         // mdpa honours mLenient rather than the narrowing options -- the same
         // "options-aware is several capabilities, not one" note as exodus. This
         // entry is what makes `--lenient` reach the C API, Fortran, Julia, R,
         // WASM and the native CLI with no per-binding code.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_mdpa(path, opts); }},
        {.mName = "medit",
         .mRead = meshioplusplus::read_medit,
         .mWrite = meshioplusplus::write_medit},
        {.mName = "medit_sol",
         .mRead = meshioplusplus::read_medit_sol,
         .mWrite = meshioplusplus::write_medit_sol},
        // A lambda: read_mfem is overloaded (grid functions). `.mesh` stays
        // medit's extension; resolve_format hands MFEM files over by content.
        {.mName = "mfem",
         .mRead = [](const std::string& path) { return meshioplusplus::read_mfem(path); },
         .mWrite = [](const std::string& path,
                      const Mesh& mesh) { meshioplusplus::write_mfem(path, mesh); },
         // MFEM honours mPiece/mPieceSet (one rank of a parallel mesh).
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_mfem(path, {}, opts);
             }},
        {.mName = "mff", .mRead = meshioplusplus::read_mff, .mWrite = meshioplusplus::write_mff},
        {.mName = "mfm",
         .mRead = meshioplusplus::read_mfm,
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_mfm(p, mm, ".16e"); }},
        {.mName = "mphbin",
         .mRead = meshioplusplus::read_mphbin,
         .mWrite = meshioplusplus::write_mphbin},
        {.mName = "mphtxt",
         .mRead = meshioplusplus::read_mphtxt,
         .mWrite = meshioplusplus::write_mphtxt},
        {.mName = "nastran",
         .mRead = meshioplusplus::read_nastran,
         .mWrite = meshioplusplus::write_nastran},
        {.mName = "nastran_op2",
         .mRead = [](const std::string& path) { return meshioplusplus::read_nastran_op2(path); },
         .mWrite = meshioplusplus::write_nastran_op2,
         // Nastran OP2 honours mTimeStep (its steps are the subcases, modes and
         // times of its result tables) and the narrowing options.
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_nastran_op2(path, opts);
             },
         .mMetadata = meshioplusplus::read_nastran_op2_metadata},
        {.mName = "netgen",
         .mRead = meshioplusplus::read_netgen,
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_netgen(p, mm, ".16e"); }},
        {.mName = "obj", .mRead = meshioplusplus::read_obj, .mWrite = meshioplusplus::write_obj},
        {.mName = "off", .mRead = meshioplusplus::read_off, .mWrite = meshioplusplus::write_off},
        // pcd/xyz take a trailing defaulted options struct, so the function pointers do
        // not convert to ReadFn -- wrapped like vtu/mdpa.
        {.mName = "pcd",
         .mRead = [](const std::string& path) { return meshioplusplus::read_pcd(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_pcd(p, mm); }},
        {.mName = "permas",
         .mRead = meshioplusplus::read_permas,
         .mWrite = meshioplusplus::write_permas},
        {.mName = "ply",
         .mRead = meshioplusplus::read_ply,
         // ply/stl default to skin=true (matching the Python shims): a volume
         // mesh writes its extracted boundary skin instead of dropping the
         // volume cells.
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_ply(p, mm, /*binary=*/true, /*skin=*/true);
             }},
        {.mName = "stl",
         .mRead = meshioplusplus::read_stl,
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_stl(p, mm, /*binary=*/false, /*skin=*/true);
             }},
        {.mName = "step",
         .mRead = [](const std::string& p) { return meshioplusplus::read_step(p); },
         // STEP is read through the optional OpenCASCADE backend (roadmap
         // §6.2.1) and written without it, one face per surface by default
         // (StepFaces::Auto: coplanar facets merged, fitted patches analytic).
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_step(p, mm); }},
        {.mName = "su2", .mRead = meshioplusplus::read_su2, .mWrite = meshioplusplus::write_su2},
        // 2-D drawings, read only as the structure meshio++ writes; both take a
        // trailing defaulted ReadOptions (strict=False is mLenient), so lambdas.
        {.mName = "gltf",
         .mRead = [](const std::string& path) { return meshioplusplus::read_gltf(path); },
         // svg/tikz are 2D-drawing formats (read only as the structure meshio++
         // writes); the flat bindings emit them with the fixed default styling
         // (per-call overrides are out of scope for v1, per registry.hpp).
         // .glb/.gltf are told apart by the writer, so the flat bindings write the
         // container the suffix names with the defaults.
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_gltf(p, mm); },
         // svg and tikz honour mLenient (strict=False: read a drawing that lacks the
         // meshio++ marker) and the narrowing options (the fill array).
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_gltf(path, GltfReadOptions{}, opts);
             }},
        // DXF R2000 drawings (roadmap 7.6.4): LINE/ARC/CIRCLE/LWPOLYLINE/SPLINE in,
        // a line mesh or a fit_curves sketch out.
        {.mName = "dxf", .mRead = meshioplusplus::read_dxf, .mWrite = meshioplusplus::write_dxf},
        {.mName = "adcirc",
         .mRead = [](const std::string& rPath) { return meshioplusplus::read_adcirc(rPath); },
         .mWrite = meshioplusplus::write_adcirc,
         .mReadEx =
             [](const std::string& rPath, const ReadOptions& rOpts) {
                 return meshioplusplus::read_adcirc(rPath, rOpts);
             },
         .mMetadata = meshioplusplus::read_adcirc_metadata},
        {.mName = "gifti",
         .mRead = meshioplusplus::read_gifti,
         .mWrite = meshioplusplus::write_gifti},
        {.mName = "freesurfer",
         .mRead = [](const std::string& rPath) { return meshioplusplus::read_freesurfer(rPath); },
         .mWrite = meshioplusplus::write_freesurfer},
        {.mName = "selafin",
         .mRead = [](const std::string& rPath) { return meshioplusplus::read_selafin(rPath); },
         .mWrite = meshioplusplus::write_selafin,
         .mReadEx =
             [](const std::string& rPath, const ReadOptions& rOpts) {
                 return meshioplusplus::read_selafin(rPath, rOpts);
             },
         .mMetadata = meshioplusplus::read_selafin_metadata},
        {.mName = "plot3d",
         .mRead = [](const std::string& rPath) { return meshioplusplus::read_plot3d(rPath); },
         .mWrite = meshioplusplus::write_plot3d},
        {.mName = "svg",
         .mRead = [](const std::string& path) { return meshioplusplus::read_svg(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_svg(p, mm); },
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_svg(path, opts); }},
        {.mName = "tikz",
         .mRead = [](const std::string& path) { return meshioplusplus::read_tikz(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_tikz(p, mm); },
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_tikz(path, opts); }},
        // A lambda, not `&read_tecplot`: the ReadOptions overload makes the
        // bare name ambiguous (the exodus/mdpa/med/cgns story again).
        {.mName = "tecplot",
         .mRead = [](const std::string& path) { return meshioplusplus::read_tecplot(path); },
         .mWrite = meshioplusplus::write_tecplot,
         // Tecplot honours mTimeStep -- selects one zone of a transient
         // (SOLUTIONTIME/STRANDID) file's timeline instead of always the
         // first. IWYU pragma: keep
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_tecplot(path, opts);
             },
         .mMetadata = meshioplusplus::read_tecplot_metadata},
        {.mName = "tetgen",
         .mRead = meshioplusplus::read_tetgen,
         .mWrite = meshioplusplus::write_tetgen},
        {.mName = "triangle",
         .mRead = meshioplusplus::read_triangle,
         .mWrite = meshioplusplus::write_triangle},
        {.mName = "ugrid",
         .mRead = meshioplusplus::read_ugrid,
         .mWrite = meshioplusplus::write_ugrid},
        {.mName = "unv",
         .mRead = [](const std::string& path) { return meshioplusplus::read_unv(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_unv(p, mm); },
         // UNV honours mTimeStep (the steps of its 2414/55/56/58 results) and the
         // narrowing options.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_unv(path, opts); },
         .mMetadata = meshioplusplus::read_unv_metadata},
        {.mName = "vti",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vti(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vti(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vti,
         .mMetadata = meshioplusplus::read_vti_metadata},
        {.mName = "vtk",
         .mRead = meshioplusplus::read_vtk,
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vtk(p, mm, /*binary=*/true, /*v51=*/true);
             }},
        {.mName = "vts",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vts(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vts(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vts,
         .mMetadata = meshioplusplus::read_vts_metadata},
        {.mName = "vtr",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtr(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vtr(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vtr,
         .mMetadata = meshioplusplus::read_vtr_metadata},
        {.mName = "vtm",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtm(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vtm(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vtm,
         .mMetadata = meshioplusplus::read_vtm_metadata},
        // The ParaView index formats take a trailing defaulted ReadOptions, so they
        // are wrapped like vtm above.
        {.mName = "pvd",
         .mRead = [](const std::string& path) { return meshioplusplus::read_pvd(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_pvd(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_pvd,
         .mMetadata = meshioplusplus::read_pvd_metadata},
        {.mName = "pvtu",
         .mRead = [](const std::string& path) { return meshioplusplus::read_pvtu(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_pvtu(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_pvtu,
         .mMetadata = meshioplusplus::read_pvtu_metadata},
        {.mName = "pvtp",
         .mRead = [](const std::string& path) { return meshioplusplus::read_pvtp(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_pvtp(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_pvtp,
         .mMetadata = meshioplusplus::read_pvtp_metadata},
        // vti/vtp/vtu take a trailing defaulted ReadOptions, so the function
        // pointers no longer convert to ReadFn -- wrapped like unv/med.
        {.mName = "vtp",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtp(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vtp(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vtp,
         .mMetadata = meshioplusplus::read_vtp_metadata},
        {.mName = "vtu",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtu(path); },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_vtu(p, mm, /*binary=*/true, kRegistryZlib);
             },
         .mReadEx = meshioplusplus::read_vtu,
         .mMetadata = meshioplusplus::read_vtu_metadata},
        {.mName = "wkt", .mRead = meshioplusplus::read_wkt, .mWrite = meshioplusplus::write_wkt},
        {.mName = "xdmf",
         .mRead = [](const std::string& path) { return meshioplusplus::read_xdmf(path); },
         // XDMF's heavy-data format follows the build: HDF companion file when
         // HDF5 is available (the Python writer's default), inline XML text
         // otherwise (the only always-available option; what WASM ships).
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
                 meshioplusplus::write_xdmf(p, mm, "HDF");
#else
                 meshioplusplus::write_xdmf(p, mm, "XML");
#endif
             },
         .mReadEx = meshioplusplus::read_xdmf,
         .mMetadata = meshioplusplus::read_xdmf_metadata},
        {.mName = "neuroglancer",
         .mRead = meshioplusplus::read_neuroglancer,
         .mWrite = meshioplusplus::write_neuroglancer},
        {.mName = "xyz",
         .mRead = [](const std::string& path) { return meshioplusplus::read_xyz(path); },
         .mWrite = [](const std::string& p, const Mesh& mm) { meshioplusplus::write_xyz(p, mm); }},
        // Side-channel info (cell-tag family names) is not carried by the flat
        // bindings -- v1 limitation, see doc/wasm.md and doc/c_api.md. Ansys
        // components travel as regions since v16.3.0, so they survive here.
        {.mName = "ansysinp",
         .mRead =
             [](const std::string& path) {
                 meshioplusplus::AnsysInfo info;
                 return meshioplusplus::read_ansysinp(path, info);
             },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::AnsysInfo info;  // no point_sets/cell_sets side channel in v1
                 meshioplusplus::write_ansysinp(p, mm, info);
             },
         // Ansys .cdb honours mLenient (skip elements with no meshio++ cell type).
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 meshioplusplus::AnsysInfo info;
                 return meshioplusplus::read_ansysinp(path, opts, info);
             }},
        {.mName = "openfoam",
         .mRead =
             [](const std::string& path) {
                 meshioplusplus::OpenFoamInfo info;
                 return meshioplusplus::read_openfoam(path, info);
             },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 // No patch-name/type side channel here, so this yields a single
                 // `defaultFaces` patch -- a valid, loadable case (see
                 // write_openfoam). The flat bindings get that rather than an error.
                 meshioplusplus::OpenFoamInfo info;
                 meshioplusplus::write_openfoam(p, mm, info);
             },
         // OpenFOAM honours mTimeStep (selects a time-directory) AND
         // mDataArrays (which fields to read) -- the OpenFoamInfo is dropped
         // here exactly as the plain reader entry drops it. IWYU pragma: keep
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 meshioplusplus::OpenFoamInfo info;
                 return meshioplusplus::read_openfoam(path, opts, info);
             },
         .mMetadata = meshioplusplus::read_openfoam_metadata},
        // Distinct from "gmsh" (4.1): 2.2 stores each element's physical tag
        // directly, so it is the version that round-trips named Cell region
        // MEMBERSHIP, not just names -- write_gmsh22 already synthesizes
        // gmsh:physical from Cell regions when the mesh has none of its own.
        // Was reachable only from Python (gmsh22_write) until this entry; the
        // flat bindings (WASM/C API/Fortran) had no way to select it at all.
        {.mName = "gmsh22",
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_gmsh22(p, mm, /*binary=*/true); }},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        // A lambda, not `&read_cgns`: the ReadOptions overload makes the bare
        // name ambiguous (the exodus/mdpa/med story again).
        {.mName = "cgns",
         .mRead = [](const std::string& path) { return meshioplusplus::read_cgns(path); },
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_cgns(p, mm, /*gzip_level=*/4); },
         // CGNS honours mTimeStep the same way exodus/med do -- selects one
         // FlowSolution_t of a transient (BaseIterativeData_t/
         // ZoneIterativeData_t) file. IWYU pragma: keep
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_cgns(path, opts); },
         .mMetadata = meshioplusplus::read_cgns_metadata},
        {.mName = "h5m",
         .mRead = meshioplusplus::read_h5m,
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 meshioplusplus::write_h5m(p, mm, /*add_global_ids=*/true, /*gzip_level=*/4);
             }},
        {.mName = "hmf",
         .mRead = meshioplusplus::read_hmf,
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_hmf(p, mm, /*gzip_level=*/4); }},
        // A lambda for the same overload reason as cgns above.
        {.mName = "vtkhdf",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtkhdf(path); },
         // UnstructuredGrid with the oldest covering version; PolyData and the
         // composite types are reached through write_vtkhdf's `Type` argument.
         .mWrite = [](const std::string& p,
                      const Mesh& mm) { meshioplusplus::write_vtkhdf(p, mm, /*gzip_level=*/4); },
         // VTKHDF honours mTimeStep (a step of a transient file), mPiece/mPieceSet
         // (one partition or composite block instead of the merged mesh), mLenient
         // (skip poly-vertex/poly-line/strip cells) and the narrowing options.
         // IWYU pragma: keep
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_vtkhdf(path, opts); },
         .mMetadata = meshioplusplus::read_vtkhdf_metadata},
        {.mName = "nastran_h5",
         .mRead = [](const std::string& path) { return meshioplusplus::read_nastran_h5(path); },
         .mWrite = meshioplusplus::write_nastran_h5,
         // MSC Nastran HDF5 honours mTimeStep (one result domain) and the
         // narrowing options, which skip the result tables not asked for.
         // IWYU pragma: keep
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_nastran_h5(path, opts);
             },
         .mMetadata = meshioplusplus::read_nastran_h5_metadata},
        {.mName = "med",
         .mRead =
             [](const std::string& path) {
                 // The family-id maps/link names/mesh metadata in MedInfo are
                 // still dropped here, but group *names* are not lost: read_med
                 // attaches them as named regions directly on the Mesh (see
                 // med_attach_point_regions/med_attach_cell_regions in med.cpp),
                 // so they reach WASM/C API/Fortran through this path too.
                 meshioplusplus::MedInfo info;
                 return meshioplusplus::read_med(path, info);
             },
         .mWrite =
             [](const std::string& p, const Mesh& mm) {
                 // No point_tags/cell_tags to hand over here (they live in the
                 // dropped MedInfo), but write_med synthesizes them from any
                 // Point/Cell regions the mesh carries (see
                 // med_point_regions_to_tags/med_cell_regions_to_tags in
                 // med.cpp), so a mesh converted from e.g. Abaqus through this
                 // registry path still carries its named groups into the file.
                 meshioplusplus::MedInfo info;
                 meshioplusplus::write_med(p, mm, info);
             },
         // MED honours `mLenient` (skip/report the enhanced `CHA` constructs
         // instead of deferring the whole file to Python) and `mTimeStep`
         // (select one step of a multi-step field) -- neither is a narrowing
         // option, the same "options-aware is several capabilities, not one"
         // note as exodus and mdpa. The `MedInfo` is dropped here exactly as
         // the plain reader entry drops it, so the flat bindings get the mesh
         // and its representable fields but not `mFieldUnits`/`mStepMeta`; that
         // is a documented gap, not a silent loss, and it is strictly more than
         // the hard failure they used to get. IWYU pragma: keep
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 meshioplusplus::MedInfo info;
                 return meshioplusplus::read_med(path, info, opts);
             },
         .mMetadata = meshioplusplus::read_med_metadata},
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        // Explicit lambda, not `&read_exodus`: the overload set now also holds
        // the ExodusInfo form, and taking its address would be ambiguous. The
        // provenance strings are dropped here, as MedInfo's are.
        {.mName = "exodus",
         .mRead = [](const std::string& path) { return meshioplusplus::read_exodus(path); },
         .mWrite = meshioplusplus::write_exodus,
         // Exodus honours mTimeStep rather than the narrowing options -- being
         // "options-aware" is not one capability but several, and a time series
         // is the one this format has. IWYU pragma: keep
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_exodus(path, opts); },
         .mMetadata = meshioplusplus::read_exodus_metadata},
        {.mName = "ugrid_nc",
         .mRead = [](const std::string& path) { return meshioplusplus::read_ugrid_nc(path); },
         .mWrite = meshioplusplus::write_ugrid_nc,
         .mReadEx =
             [](const std::string& path, const ReadOptions& opts) {
                 return meshioplusplus::read_ugrid_nc(path, opts);
             },
         .mMetadata = meshioplusplus::read_ugrid_nc_metadata},
#endif
#ifdef MESHIOPLUSPLUS_HAS_ADIOS2
        {.mName = "vtx",
         .mRead = [](const std::string& path) { return meshioplusplus::read_vtx(path); },
         .mWrite = meshioplusplus::write_vtx,
         // VTX honours mTimeStep, mPointsOnly, mDataArrays and mGhosts.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_vtx(path, opts); },
         .mMetadata = meshioplusplus::read_vtx_metadata},
#endif
#ifdef MESHIOPLUSPLUS_HAS_TECIO
        {.mName = "szplt",
         .mRead = [](const std::string& path) { return meshioplusplus::read_szplt(path); },
         .mWrite = meshioplusplus::write_szplt,
         // .szplt honours mTimeStep (the distinct solution times), as tecplot.
         .mReadEx = [](const std::string& path,
                       const ReadOptions& opts) { return meshioplusplus::read_szplt(path, opts); },
         .mMetadata = meshioplusplus::read_szplt_metadata},
#endif
    };
    return rows;
}

// The rows' entries of one kind, as the map type its accessor returns.
template <class TMap, class TField>
TMap registry_collect(TField RegistryRow::*pField) {
    TMap out;
    for (const RegistryRow& r_row : registry_rows())
        if (r_row.*pField)
            out.emplace(r_row.mName, r_row.*pField);
    return out;
}

}  // namespace

const std::map<std::string, ReadFn>& registry_readers() {
    static const std::map<std::string, ReadFn> m = registry_guard_readers(
        registry_collect<std::map<std::string, ReadFn>>(&RegistryRow::mRead));
    return m;
}

const std::map<std::string, WriteFn>& registry_writers() {
    static const std::map<std::string, WriteFn> m =
        registry_collect<std::map<std::string, WriteFn>>(&RegistryRow::mWrite);
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
        {".fil", "abaqus_fil"},
        {".cdb", "ansysinp"},
        {".frd", "frd"},
        {".k", "lsdyna"},
        {".key", "lsdyna"},
        {".dyn", "lsdyna"},
        {".mail", "code_aster"},
        {".pat", "patran"},
        {".out", "patran"},
        {".neu", "femap"},
        {".xda", "libmesh"},
        {".xdr", "libmesh"},
        {".xda.gz", "libmesh"},
        {".xdr.gz", "libmesh"},
        {".xda.bz2", "libmesh"},
        {".xdr.bz2", "libmesh"},
        {".rad", "radioss"},
        {".t19", "marc_t19"},
        {".feb", "febio"},
        {".xplt", "xplt"},
        {".rst", "ansys_rst"},
        {".rth", "ansys_rst"},
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
        {".meshb", "medit"},
        {".sol", "medit_sol"},
        {".solb", "medit_sol"},
        {".mfm", "mfm"},
        {".mphbin", "mphbin"},
        {".mphtxt", "mphtxt"},
        {".bdf", "nastran"},
        {".nas", "nastran"},
        {".fem", "nastran"},
        {".op2", "nastran_op2"},
        {".vol", "netgen"},
        {".vol.gz", "netgen"},
        {".obj", "obj"},
        // OpenFOAM: the `.foam` marker file. A case *directory* has no
        // extension at all, so that form still needs an explicit format.
        {".foam", "openfoam"},
        {".off", "off"},
        {".pcd", "pcd"},
        {".post", "permas"},
        {".dato", "permas"},
        // GiD postprocess. All four are compound extensions, and
        // resolve_format() tries the longest suffix first, so ".post.msh"
        // resolves to "gid" rather than falling through to ".msh" -> gmsh.
        // Registered even when the build has no gidpost: write_gid() then
        // throws naming the missing build flags, which is strictly better
        // than ".post.msh" silently resolving to another format.
        {".post.msh", "gid"},
        // Kratos CAD JSON; a bare ".json" stays unclaimed (sniffed by content).
        {".cad.json", "cad_json"},
        {".post.res", "gid"},
        {".post.bin", "gid"},
        {".post.h5", "gid"},
        {".ply", "ply"},
        {".igs", "iges"},
        {".iges", "iges"},
        {".step", "step"},
        {".stl", "stl"},
        {".stp", "step"},
        {".su2", "su2"},
        {".glb", "gltf"},
        {".gltf", "gltf"},
        {".dxf", "dxf"},
        {".14", "adcirc"},
        {".gii", "gifti"},
        {".pial", "freesurfer"},
        {".white", "freesurfer"},
        {".inflated", "freesurfer"},
        {".orig", "freesurfer"},
        {".smoothwm", "freesurfer"},
        {".sphere", "freesurfer"},
        {".slf", "selafin"},
        {".srf", "selafin"},
        {".grd", "adcirc"},
        {".p3d", "plot3d"},
        {".p3dfmt", "plot3d"},
        {".x", "plot3d"},
        {".g", "plot3d"},
        {".svg", "svg"},
        {".tikz", "tikz"},
        // .node/.ele stay with tetgen for backward compatibility (3D pairs
        // are the common case for the flat bindings); reading Triangle 2D
        // .node/.ele files there needs an explicit format="triangle". Only
        // .poly defaults to triangle.
        {".dat", "tecplot"},
        {".tec", "tecplot"},
        {".plt", "tecplot"},
        {".ele", "tetgen"},
        {".node", "tetgen"},
        {".smesh", "tetgen"},
        {".poly", "triangle"},
        {".ugrid", "ugrid"},
        {".unv", "unv"},
        {".uff", "unv"},
        {".vti", "vti"},
        {".vtk", "vtk"},
        {".vts", "vts"},
        {".vtr", "vtr"},
        {".vtm", "vtm"},
        {".pvd", "pvd"},
        {".pvtu", "pvtu"},
        {".pvtp", "pvtp"},
        {".vtp", "vtp"},
        {".vtu", "vtu"},
        {".wkt", "wkt"},
        {".xdmf", "xdmf"},
        {".xmf", "xdmf"},
        // Point-cloud text: a convention rather than a specification, so every alias the
        // wild uses maps to the one column-sniffing reader.
        {".xyz", "xyz"},
        {".xyzn", "xyz"},
        {".xyzrgb", "xyz"},
        {".asc", "xyz"},
        {".pts", "xyz"},
        {".txt", "xyz"},
        {".msh", "gmsh"},
        {".cgns", "cgns"},
        {".h5m", "h5m"},
        {".hmf", "hmf"},
        {".vtkhdf", "vtkhdf"},
        {".hdf", "vtkhdf"},
        {".med", "med"},
        {".h5", "nastran_h5"},
        {".bp", "vtx"},
        {".szplt", "szplt"},
        {".e", "exodus"},
        {".exo", "exodus"},
        {".ex2", "exodus"},
        {".nc", "ugrid_nc"},
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
// Whether `rPath` exists and opens with an MFEM mesh header (`MFEM mesh v1.x`,
// `MFEM NC mesh ...`); only for the `.mesh` suffix Medit and MFEM share.
bool registry_is_mfem_mesh(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        return false;
    char buf[64];
    in.read(buf, sizeof(buf));
    std::string head(buf, static_cast<std::size_t>(in.gcount()));
    const std::size_t first = head.find_first_not_of(" \t\r\n");
    return first != std::string::npos && head.compare(first, 5, "MFEM ") == 0;
}

// Whether `rPath` exists and opens as a Marc input deck; only for the `.dat`
// suffix Tecplot and Marc share.
bool registry_is_marc_dat(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        return false;
    std::string head(65536, '\0');
    in.read(head.data(), static_cast<std::streamsize>(head.size()));
    head.resize(static_cast<std::size_t>(in.gcount()));
    return is_marc_deck(head);
}

// `d3plot01`, `d3plot02`... beside an existing `d3plot` (or `d3part01`... beside a `d3part`).
bool registry_is_d3plot_member(const std::string& rPath) {
    namespace fs = std::filesystem;
    const fs::path path(rPath);
    std::string name = path.filename().string();
    for (char& c : name)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (name.size() <= 6 ||
        (name.compare(0, 6, "d3plot") != 0 && name.compare(0, 6, "d3part") != 0))
        return false;
    for (std::size_t i = 6; i < name.size(); ++i)
        if (name[i] < '0' || name[i] > '9')
            return false;
    std::error_code ec;
    return fs::is_regular_file(path.parent_path() / path.filename().string().substr(0, 6), ec);
}

std::string basename_of(const std::string& rPath) {
    auto pos = rPath.find_last_of("/\\");
    return pos == std::string::npos ? rPath : rPath.substr(pos + 1);
}

/// `resolve_format`, looking at an existing file's content only when
/// `Content` is set: a write never lets the file it replaces choose its format.
std::string registry_resolve(const std::string& rPath, const std::string& rFormat, bool Content) {
    if (!rFormat.empty())
        return rFormat;
    const auto& defaults = registry_extension_defaults();
    const std::string base = basename_of(rPath);
    // Z88's files have fixed names; `.txt` alone is xyz's.
    if (is_z88_filename(base))
        return "z88";
    // LS-DYNA's state database has a fixed name too, and no extension; its
    // numbered members (`d3plot01`...) go to the reader, which names the base.
    if (is_d3plot_filename(base) || registry_is_d3plot_member(rPath))
        return "lsdyna_d3plot";
    if (is_binout_filename(base))
        return "lsdyna_binout";
    // OpenRadioss animation files: `<run>A001`..., no extension.
    if (base.find('.') == std::string::npos && is_radioss_anim_filename(base))
        return "radioss_anim";
    // ... and their time-history files: `<run>T01`...
    if (is_radioss_th_filename(base))
        return "radioss_th";
    for (std::size_t pos = base.find('.'); pos != std::string::npos;
         pos = base.find('.', pos + 1)) {
        const std::string suffix = base.substr(pos);
        auto it = defaults.find(suffix);
        if (it == defaults.end())
            continue;
        // `.mesh` is both Medit's and MFEM's: an existing file whose first line
        // names an MFEM mesh goes to mfem.
        if (Content && suffix == ".mesh" && registry_is_mfem_mesh(rPath))
            return "mfem";
        // `.dat` is Tecplot's, and Marc's input deck's when it opens as one.
        if (Content && suffix == ".dat" && registry_is_marc_dat(rPath))
            return "marc";
        return it->second;
    }
    // `cube.msh.gz`: the format is the one the name had before compression.
    const detail::Compression kind = registry_compression(rPath);
    if (kind != detail::Compression::None)
        return registry_resolve(detail::decompressed_name(rPath, kind), "", /*Content=*/false);
    throw meshioplusplus::ReadError("meshio++: cannot infer format from '" + rPath +
                                    "' -- pass an explicit format argument");
}

}  // namespace

std::string resolve_format(const std::string& rPath, const std::string& rFormat) {
    return registry_resolve(rPath, rFormat, /*Content=*/true);
}

std::string resolve_write_format(const std::string& rPath, const std::string& rFormat) {
    return registry_resolve(rPath, rFormat, /*Content=*/false);
}

const std::unordered_map<std::string, ReadExFn>& registry_readers_ex() {
    // Sparse by design -- populated per format as native selective-read support
    // lands. An absent format falls back to a full read in registry_read().
    static const std::unordered_map<std::string, ReadExFn> m =
        registry_collect<std::unordered_map<std::string, ReadExFn>>(&RegistryRow::mReadEx);
    return m;
}

const std::unordered_map<std::string, MetadataFn>& registry_metadata_readers() {
    static const std::unordered_map<std::string, MetadataFn> m =
        registry_collect<std::unordered_map<std::string, MetadataFn>>(&RegistryRow::mMetadata);
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
    const detail::Compression kind = registry_compression(rPath);
    if (kind != detail::Compression::None) {
        registry_refuse_multi_file(rFormat, rPath);
        const detail::DecompressedFile plain(rPath, kind);
        return registry_read(plain.Path(), rFormat, rOptions);
    }
    auto it = registry_readers_ex().find(rFormat);
    if (it != registry_readers_ex().end())
        return detail::guarded_read(rFormat.c_str(), [&] { return it->second(rPath, rOptions); });
    // No native selective path: a full read is still the correct answer.
    return registry_full_reader(rFormat)(rPath);
}

MeshMetadata registry_read_metadata(const std::string& rPath, const std::string& rFormat,
                                    const ReadOptions& rOptions) {
    const detail::Compression kind = registry_compression(rPath);
    if (kind != detail::Compression::None) {
        registry_refuse_multi_file(rFormat, rPath);
        const detail::DecompressedFile plain(rPath, kind);
        return registry_read_metadata(plain.Path(), rFormat, rOptions);
    }
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
            MeshMetadata meta =
                detail::guarded_read(rFormat.c_str(), [&] { return it->second(rPath, rOptions); });
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
    if (rFormat == "cgns" || rFormat == "h5m" || rFormat == "hmf" || rFormat == "med" ||
        rFormat == "vtkhdf" || rFormat == "nastran_h5")
        return "HDF5";
#endif
#ifndef MESHIOPLUSPLUS_HAS_NETCDF
    if (rFormat == "exodus" || rFormat == "ugrid_nc")
        return "netCDF";
#endif
#ifndef MESHIOPLUSPLUS_HAS_ADIOS2
    if (rFormat == "vtx")
        return "ADIOS2";
#endif
#ifndef MESHIOPLUSPLUS_HAS_TECIO
    if (rFormat == "szplt")
        return "TecIO";
#endif
    // Deliberately NOT an arm for cgns/cgnslib: the format is fully readable
    // and writable without it (the hand-rolled ADF-over-HDF5 path), so
    // reporting it "compiled out" would be wrong. The cgnslib-only
    // capabilities -- ADF containers and the CGNS 3.x NGON_n/NFACE_n layout --
    // report themselves by name from read_cgns_mll instead.
    return nullptr;
}

}  // namespace meshioplusplus
