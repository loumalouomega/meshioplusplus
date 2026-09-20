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
// VTK XML parallel indices `.pvtu` (PUnstructuredGrid) and `.pvtp` (PPolyData).
// The two differ only in the root tag and the piece format, so one set of
// `pvtu_`-prefixed helpers serves both and the public `_pvtp` functions in
// `formats/pvtp.hpp` are defined at the bottom of this file.
//
// Writing means validating that every piece declares identical arrays, deriving
// `vtkGhostType` from `partition:ghost`, then delegating each piece to
// write_vtu_codec / write_vtp_codec; reading means parsing the index and
// delegating every piece to read_vtu / read_vtp (see formats/pindex_common.hpp).
//
// Anonymous-namespace helpers are prefixed `pvtu_`, as the amalgamation
// concatenates every translation unit.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "pindex_common.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/pvtp.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/crop.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

enum class pvtu_kind { Pvtu, Pvtp };

struct pvtu_kind_info {
    const char* mFormat;
    const char* mRoot;
    const char* mExt;
};

const pvtu_kind_info& pvtu_info(pvtu_kind Kind) {
    static const pvtu_kind_info sPvtu{"pvtu", "PUnstructuredGrid", ".vtu"};
    static const pvtu_kind_info sPvtp{"pvtp", "PPolyData", ".vtp"};
    return Kind == pvtu_kind::Pvtu ? sPvtu : sPvtp;
}

// vtkDataSetAttributes: DUPLICATECELL = 1, DUPLICATEPOINT = 1. REFINEDCELL (8) and
// HIDDENCELL/HIDDENPOINT (32) are read back faithfully; nothing here produces them.
constexpr std::uint8_t pvtu_duplicate_cell = 1;
constexpr std::uint8_t pvtu_duplicate_point = 1;

// --- declarations ------------------------------------------------------------

struct pvtu_decl {
    std::string mType;
    std::size_t mComp = 1;

    bool operator==(const pvtu_decl& rOther) const {
        return mType == rOther.mType && mComp == rOther.mComp;
    }
    bool operator!=(const pvtu_decl& rOther) const { return !(*this == rOther); }
};

struct pvtu_decls {
    pvtu_decl mPoints;
    std::map<std::string, pvtu_decl> mPointData;
    std::map<std::string, pvtu_decl> mCellData;
};

pvtu_decl pvtu_decl_of(const NDArray& rArray) {
    pvtu_decl d;
    d.mType = detail::vtu_type_str(rArray.Dtype());
    d.mComp = rArray.Shape().size() == 2 ? rArray.Shape()[1] : 1;
    return d;
}

std::string pvtu_fmt_decl(const pvtu_decl& rDecl) {
    return rDecl.mType + " with " + std::to_string(rDecl.mComp) +
           (rDecl.mComp == 1 ? " component" : " components");
}

std::size_t pvtu_num_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += cb.NumCells();
    return n;
}

// What a piece will hold, as `<PDataArray>` declares it.
pvtu_decls pvtu_declarations_of(const Mesh& rMesh) {
    pvtu_decls d;
    d.mPoints.mType = detail::vtu_type_str(rMesh.Points().Dtype());
    d.mPoints.mComp = 3;
    for (const std::string& name : rMesh.PointDataNames())
        d.mPointData[name] = pvtu_decl_of(rMesh.PointData(name));
    for (const std::string& name : rMesh.CellDataNames()) {
        bool have = false;
        pvtu_decl decl;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b) {
            const pvtu_decl cur = pvtu_decl_of(rMesh.CellData(name, b));
            if (!have) {
                decl = cur;
                have = true;
            } else if (cur != decl) {
                throw WriteError("meshio++: cell_data '" + name +
                                 "' has differing types or component counts across cell blocks");
            }
        }
        if (have)
            d.mCellData[name] = decl;
    }
    return d;
}

// Refuse a piece list whose pieces do not declare identical arrays, before
// anything is emitted. Piece 0 is the reference. A piece with no cells is exempt
// from the cell_data comparison only when it carries no cell arrays (an idle rank
// legitimately allocated nothing); points are never exempt.
pvtu_decls pvtu_validate(pvtu_kind Kind, const std::vector<const Mesh*>& rPieces) {
    const std::string format = pvtu_info(Kind).mFormat;
    std::vector<pvtu_decls> decls;
    decls.reserve(rPieces.size());
    for (const Mesh* p : rPieces)
        decls.push_back(pvtu_declarations_of(*p));

    const std::string tail =
        "; every piece of a parallel index must declare identical arrays (name, type, "
        "NumberOfComponents)";
    for (std::size_t k = 1; k < decls.size(); ++k) {
        const std::string head = "meshio++: " + format + ": piece " + std::to_string(k);
        const pvtu_decls& ref = decls[0];
        const pvtu_decls& cur = decls[k];
        if (cur.mPoints != ref.mPoints)
            throw WriteError(head + " declares Points as " + pvtu_fmt_decl(cur.mPoints) +
                             ", but piece 0 declares it as " + pvtu_fmt_decl(ref.mPoints) + tail);

        struct section {
            const char* mLabel;
            const std::map<std::string, pvtu_decl>* mRef;
            const std::map<std::string, pvtu_decl>* mCur;
        };
        const section sections[] = {{"point_data", &ref.mPointData, &cur.mPointData},
                                    {"cell_data", &ref.mCellData, &cur.mCellData}};
        for (const section& s : sections) {
            if (std::string(s.mLabel) == "cell_data" && s.mCur->empty() &&
                pvtu_num_cells(*rPieces[k]) == 0)
                continue;
            std::set<std::string> names;
            for (const auto& kv : *s.mRef)
                names.insert(kv.first);
            for (const auto& kv : *s.mCur)
                names.insert(kv.first);
            for (const std::string& name : names) {
                const auto in_cur = s.mCur->find(name);
                const auto in_ref = s.mRef->find(name);
                if (in_cur == s.mCur->end())
                    throw WriteError(head + " is missing " + s.mLabel + " '" + name +
                                     "', which piece 0 declares" + tail);
                if (in_ref == s.mRef->end())
                    throw WriteError(head + " declares " + s.mLabel + " '" + name +
                                     "', which piece 0 does not" + tail);
                if (in_ref->second != in_cur->second)
                    throw WriteError(head + " declares " + s.mLabel + " '" + name + "' as " +
                                     pvtu_fmt_decl(in_cur->second) +
                                     ", but piece 0 declares it as " +
                                     pvtu_fmt_decl(in_ref->second) + tail);
            }
        }
    }
    return decls.front();
}

// --- ghosts ------------------------------------------------------------------

std::int64_t pvtu_ghost_level(const std::vector<const Mesh*>& rPieces) {
    std::int64_t level = 0;
    for (const Mesh* p : rPieces) {
        if (!p->HasCellData(pidx::kPartitionGhost))
            continue;
        for (std::size_t b = 0; b < p->CellDataNumBlocks(pidx::kPartitionGhost); ++b)
            for (const std::int64_t v : detail::vtu_to_int64(p->CellData(pidx::kPartitionGhost, b)))
                level = std::max(level, v);
    }
    return level;
}

bool pvtu_needs_ghosts(const Mesh& rMesh) {
    return rMesh.HasCellData(pidx::kPartitionGhost) &&
           !(rMesh.HasCellData(pidx::kGhostName) && rMesh.HasPointData(pidx::kGhostName));
}

// A copy of @p rMesh carrying `vtkGhostType` derived from `partition:ghost`.
// Cells: 0 for an owned cell (layer 0), DUPLICATECELL for any halo layer. Points:
// DUPLICATEPOINT when no owned cell of this piece references the point. An array
// the caller already supplied is passed through unchanged.
Mesh pvtu_derive_ghosts(const Mesh& rMesh) {
    Mesh out = detail::clone_mesh(rMesh);
    const bool have_cell = rMesh.HasCellData(pidx::kGhostName);
    const bool have_point = rMesh.HasPointData(pidx::kGhostName);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<std::uint8_t> owned(rMesh.NumPoints(), 0);
    std::vector<NDArray> cell_gh;
    cell_gh.reserve(nblocks);

    for (std::size_t b = 0; b < nblocks; ++b) {
        const std::vector<std::int64_t> layers =
            detail::vtu_to_int64(rMesh.CellData(pidx::kPartitionGhost, b));
        const auto cb = rMesh.Cells(b);
        if (layers.size() != cb.NumCells())
            throw WriteError("meshio++: partition:ghost has " + std::to_string(layers.size()) +
                             " values for a block of " + std::to_string(cb.NumCells()) + " cells");
        NDArray gh = NDArray::Uninit(DType::UInt8, {layers.size()});
        std::uint8_t* g = gh.As<std::uint8_t>();
        for (std::size_t c = 0; c < layers.size(); ++c)
            g[c] = layers[c] > 0 ? pvtu_duplicate_cell : 0;
        cell_gh.push_back(std::move(gh));

        auto own = [&](std::int64_t Node) {
            if (Node >= 0 && static_cast<std::size_t>(Node) < owned.size())
                owned[static_cast<std::size_t>(Node)] = 1;
        };
        if (cb.IsPolyhedron()) {
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                if (layers[c] != 0)
                    continue;
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    for (std::size_t i = 0; i < face.second; ++i)
                        own(face.first[i]);
                }
            }
        } else if (cb.IsRagged()) {
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                if (layers[c] == 0)
                    for (std::size_t i = 0; i < cb.RowSize(c); ++i)
                        own(cb.Row(c)[i]);
        } else {
            const std::vector<std::int64_t> conn = detail::vtu_to_int64(cb.Conn());
            const std::size_t npc = cb.NodesPerCell();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                if (layers[c] == 0)
                    for (std::size_t i = 0; i < npc; ++i)
                        own(conn[c * npc + i]);
        }
    }

    if (!have_cell)
        out.AddCellData(pidx::kGhostName, std::move(cell_gh));
    if (!have_point) {
        NDArray pg = NDArray::Uninit(DType::UInt8, {owned.size()});
        std::uint8_t* g = pg.As<std::uint8_t>();
        for (std::size_t i = 0; i < owned.size(); ++i)
            g[i] = owned[i] ? 0 : pvtu_duplicate_point;
        out.AddPointData(pidx::kGhostName, std::move(pg));
    }
    return out;
}

// --- carving a single mesh into parts ----------------------------------------

// One piece per part id 0..max of the integer cell_data @p rKey. Parts with no
// cells are kept (they are real ranks; skipping one would renumber the rest).
// Each piece keeps every block (possibly empty) and only the points its cells use.
std::vector<Mesh> pvtu_carve(const Mesh& rMesh, const std::string& rKey) {
    std::int64_t max_id = -1;
    for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(rKey); ++b) {
        const NDArray& ids = rMesh.CellData(rKey, b);
        if (ids.Dtype() == DType::Float32 || ids.Dtype() == DType::Float64)
            throw WriteError("meshio++: cell_data '" + rKey + "' must hold integer part ids");
        for (const std::int64_t v : detail::vtu_to_int64(ids)) {
            if (v < 0)
                throw WriteError("meshio++: cell_data '" + rKey + "' holds a negative part id");
            max_id = std::max(max_id, v);
        }
    }
    const std::int64_t nparts = std::max<std::int64_t>(max_id + 1, 1);
    std::vector<Mesh> pieces;
    pieces.reserve(static_cast<std::size_t>(nparts));
    for (std::int64_t p = 0; p < nparts; ++p)
        pieces.push_back(std::move(
            crop_predicate(rMesh, rKey, RefineCompare::Equal, static_cast<double>(p)).mMesh));
    return pieces;
}

// --- write -------------------------------------------------------------------

void pvtu_write_pieces(pvtu_kind Kind, const std::string& rPath,
                       const std::vector<const Mesh*>& rPieces, bool Binary,
                       detail::VtkCodec Codec) {
    const pvtu_kind_info& info = pvtu_info(Kind);
    if (rPieces.empty())
        throw WriteError(std::string("meshio++: ") + info.mFormat +
                         ": a parallel index needs at least one piece");
    if (Binary && Codec != detail::VtkCodec::None)
        detail::vtk_codec_require_write(Codec);

    // Ghost arrays are derived on a copy: the caller's meshes are never mutated.
    std::vector<Mesh> derived;
    derived.reserve(rPieces.size());
    std::vector<const Mesh*> prepared;
    prepared.reserve(rPieces.size());
    for (const Mesh* p : rPieces) {
        if (pvtu_needs_ghosts(*p)) {
            derived.push_back(pvtu_derive_ghosts(*p));
            prepared.push_back(&derived.back());
        } else {
            prepared.push_back(p);
        }
    }

    // Validate before create_directories, so a refusal leaves no half-written tree.
    const pvtu_decls decls = pvtu_validate(Kind, prepared);
    const std::int64_t level = pvtu_ghost_level(rPieces);

    const fs::path index_path(rPath);
    const std::string stem = index_path.stem().string();
    const fs::path dir =
        index_path.parent_path().empty() ? fs::path(stem) : index_path.parent_path() / stem;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec)
        throw WriteError("Could not create directory for " + std::string(info.mExt) +
                         " pieces: " + dir.string() + ": " + ec.message());

    const std::vector<std::string> names = pidx::piece_names(stem, prepared.size(), info.mExt);
    const detail::VtkCodec piece_codec = Binary ? Codec : detail::VtkCodec::None;
    for (std::size_t i = 0; i < prepared.size(); ++i) {
        const std::string path = (dir / names[i]).string();
        if (Kind == pvtu_kind::Pvtu)
            write_vtu_codec(path, *prepared[i], Binary, piece_codec);
        else
            write_vtp_codec(path, *prepared[i], Binary, piece_codec);
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"" << info.mRoot << "\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<" << info.mRoot << " GhostLevel=\"" << level << "\">\n";
    const std::pair<const char*, const std::map<std::string, pvtu_decl>*> sections[] = {
        {"PPointData", &decls.mPointData}, {"PCellData", &decls.mCellData}};
    for (const auto& section : sections) {
        if (section.second->empty())
            continue;
        os << "<" << section.first << ">\n";
        for (const auto& kv : *section.second) {
            os << "<PDataArray type=\"" << kv.second.mType << "\" Name=\""
               << pidx::escape_attr(kv.first) << "\"";
            if (kv.second.mComp != 1)
                os << " NumberOfComponents=\"" << kv.second.mComp << "\"";
            os << "/>\n";
        }
        os << "</" << section.first << ">\n";
    }
    os << "<PPoints>\n<PDataArray type=\"" << decls.mPoints.mType
       << "\" Name=\"Points\" NumberOfComponents=\"3\"/>\n</PPoints>\n";
    for (const std::string& name : names)
        os << "<Piece Source=\"" << pidx::escape_attr(stem + "/" + name) << "\"/>\n";
    os << "</" << info.mRoot << ">\n";
    os << "</VTKFile>\n";
    if (!os)
        throw WriteError("Failed while writing: " + rPath);
}

void pvtu_write_mesh(pvtu_kind Kind, const std::string& rPath, const Mesh& rMesh, bool Binary,
                     detail::VtkCodec Codec, const std::string& rPartKey) {
    if (rPartKey.empty() || !rMesh.HasCellData(rPartKey)) {
        pvtu_write_pieces(Kind, rPath, {&rMesh}, Binary, Codec);
        return;
    }
    const std::vector<Mesh> parts = pvtu_carve(rMesh, rPartKey);
    std::vector<const Mesh*> ptrs;
    ptrs.reserve(parts.size());
    for (const Mesh& m : parts)
        ptrs.push_back(&m);
    pvtu_write_pieces(Kind, rPath, ptrs, Binary, Codec);
}

// --- read --------------------------------------------------------------------

std::vector<fs::path> pvtu_parse_index(pvtu_kind Kind, const std::string& rPath) {
    const pvtu_kind_info& info = pvtu_info(Kind);
    const std::string head = std::string("meshio++: ") + info.mFormat + ": ";
    pugi::xml_document doc;
    const pugi::xml_parse_result parsed = doc.load_file(rPath.c_str());
    if (!parsed)
        throw ReadError(head + "could not parse " + rPath + ": " + parsed.description());
    const pugi::xml_node root = doc.child("VTKFile");
    if (!root)
        throw ReadError(head + "expected tag 'VTKFile': " + rPath);
    const std::string type = root.attribute("type").as_string();
    if (type != info.mRoot)
        throw ReadError(head + "expected type " + info.mRoot + ", got '" + type + "': " + rPath);
    const pugi::xml_node body = root.child(info.mRoot);
    if (!body)
        throw ReadError(head + "expected tag '" + info.mRoot + "': " + rPath);

    std::vector<fs::path> sources;
    for (const pugi::xml_node piece : body.children("Piece"))
        sources.push_back(pidx::resolve_path(rPath, piece.attribute("Source").as_string(""),
                                             info.mFormat, "Source"));
    return sources;
}

Mesh pvtu_read_one(const fs::path& rSource, const ReadOptions& rOpts, GhostPolicy Ghosts) {
    Mesh piece = pidx::read_child(rSource, rOpts, PvtuReadOptions{Ghosts}, /*Wide=*/false);
    if (Ghosts == GhostPolicy::Drop)
        return pidx::drop_ghosts(std::move(piece));
    return piece;
}

Mesh pvtu_read(pvtu_kind Kind, const std::string& rPath, const ReadOptions& rOpts,
               const PvtuReadOptions& rGhost) {
    const std::vector<fs::path> sources = pvtu_parse_index(Kind, rPath);
    if (sources.empty())
        return pidx::empty_mesh();
    if (rOpts.mPieceSet)
        return pvtu_read_one(sources[rOpts.ResolvePiece(sources.size())], rOpts, rGhost.mGhosts);

    std::vector<Mesh> pieces;
    std::vector<std::string> names;
    pieces.reserve(sources.size());
    for (std::size_t i = 0; i < sources.size(); ++i) {
        pieces.push_back(pvtu_read_one(sources[i], rOpts, rGhost.mGhosts));
        names.push_back("piece_" + std::to_string(i));
    }
    return pidx::merge_pieces(std::move(pieces), names);
}

MeshMetadata pvtu_read_metadata(pvtu_kind Kind, const std::string& rPath,
                                const ReadOptions& rOpts) {
    const std::vector<fs::path> sources = pvtu_parse_index(Kind, rPath);
    std::vector<MeshMetadata> parts;
    parts.reserve(sources.size());
    for (const fs::path& s : sources)
        parts.push_back(pidx::read_child_metadata(s, rOpts, /*Wide=*/false));
    return pidx::aggregate_metadata(parts, pvtu_info(Kind).mFormat);
}

detail::VtkCodec pvtu_zlib_codec(bool Zlib) {
    return Zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None;
}

}  // namespace

// --- .pvtu -------------------------------------------------------------------

void write_pvtu(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_pvtu_codec(rPath, rMesh, binary, pvtu_zlib_codec(zlib));
}

void write_pvtu_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                      detail::VtkCodec codec, const std::string& rPartKey) {
    pvtu_write_mesh(pvtu_kind::Pvtu, rPath, rMesh, binary, codec, rPartKey);
}

void write_pvtu_pieces_codec(const std::string& rPath, const std::vector<const Mesh*>& rPieces,
                             bool binary, detail::VtkCodec codec) {
    pvtu_write_pieces(pvtu_kind::Pvtu, rPath, rPieces, binary, codec);
}

Mesh read_pvtu(const std::string& rPath, const ReadOptions& rOpts, const PvtuReadOptions& rGhost) {
    return pvtu_read(pvtu_kind::Pvtu, rPath, rOpts, rGhost);
}

MeshMetadata read_pvtu_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    return pvtu_read_metadata(pvtu_kind::Pvtu, rPath, rOpts);
}

// --- .pvtp -------------------------------------------------------------------

void write_pvtp(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_pvtp_codec(rPath, rMesh, binary, pvtu_zlib_codec(zlib));
}

void write_pvtp_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                      detail::VtkCodec codec, const std::string& rPartKey) {
    pvtu_write_mesh(pvtu_kind::Pvtp, rPath, rMesh, binary, codec, rPartKey);
}

void write_pvtp_pieces_codec(const std::string& rPath, const std::vector<const Mesh*>& rPieces,
                             bool binary, detail::VtkCodec codec) {
    pvtu_write_pieces(pvtu_kind::Pvtp, rPath, rPieces, binary, codec);
}

Mesh read_pvtp(const std::string& rPath, const ReadOptions& rOpts, const PvtuReadOptions& rGhost) {
    return pvtu_read(pvtu_kind::Pvtp, rPath, rOpts, rGhost);
}

MeshMetadata read_pvtp_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    return pvtu_read_metadata(pvtu_kind::Pvtp, rPath, rOpts);
}

}  // namespace meshioplusplus
