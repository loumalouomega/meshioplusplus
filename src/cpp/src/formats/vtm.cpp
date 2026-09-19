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
// VTK XML MultiBlock. The index carries no geometry at all -- writing means
// carving the mesh into one piece per CellBlock and delegating every piece to
// write_vtu_codec; reading means parsing the index and delegating every piece
// to read_vtu/read_vtp, then combining with operations/merge.hpp.
//
// Anonymous-namespace helpers are prefixed `vtm_`, as the amalgamation
// concatenates every translation unit and `vtu_`/`vts_`/`vtr_` are taken.

// System includes
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/vtm.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

// Build a standalone single-block mesh: block `Idx`'s own cells (whichever
// representation -- rectangular, polygon, or polyhedron), the mesh's full
// point array (pruned below), and that block's slice of every cell_data key.
// field_data is not carried onto pieces: it describes the whole mesh, not one
// block, and .vtm has no per-block place to put it (see vtm.md).
Mesh vtm_extract_piece(const Mesh& rMesh, std::size_t Idx) {
    Mesh piece;
    piece.AssignPoints(detail::data_owned_copy(rMesh.Points()));

    const auto cb = rMesh.Cells(Idx);
    if (cb.IsPolyhedron()) {
        std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            cells[r].resize(cb.NumFaces(r));
            for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                const auto face = cb.Face(r, f);
                cells[r][f].assign(face.first, face.first + face.second);
            }
        }
        piece.AddPolyhedronBlock(std::string(cb.Type()), std::move(cells));
    } else if (cb.IsRagged()) {
        std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
        for (std::size_t r = 0; r < cb.NumCells(); ++r)
            rows[r].assign(cb.Row(r), cb.Row(r) + cb.RowSize(r));
        piece.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
    } else {
        piece.AddCellBlock(std::string(cb.Type()), detail::data_owned_copy(cb.Conn()));
    }

    for (const std::string& name : rMesh.PointDataNames())
        piece.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.CellDataNames())
        if (rMesh.CellDataNumBlocks(name) > Idx)
            piece.AddCellData(name, {detail::data_owned_copy(rMesh.CellData(name, Idx))});

    return piece;
}

// Every `<DataSet>` under `pNode`, in document order, regardless of how many
// `<Block>` levels separate them from it (xpath is not linked in, so this is
// a plain recursive walk rather than a `"//DataSet"` select).
void vtm_collect_datasets(const pugi::xml_node& rNode, std::vector<pugi::xml_node>& rOut) {
    for (pugi::xml_node child : rNode.children()) {
        if (std::string(child.name()) == "DataSet")
            rOut.push_back(child);
        else
            vtm_collect_datasets(child, rOut);
    }
}

// One piece, resolved from its `<DataSet>` element: the path to read (relative
// to the index file's own directory) and the name to give it.
struct vtm_piece_ref {
    fs::path mPath;
    std::string mName;
};

std::vector<vtm_piece_ref> vtm_parse_index(const std::string& rPath, pugi::xml_document& rDoc) {
    const pugi::xml_parse_result parsed = rDoc.load_file(rPath.c_str());
    if (!parsed)
        throw ReadError("Could not parse .vtm XML: " + rPath + ": " + parsed.description());

    const pugi::xml_node root = rDoc.child("VTKFile");
    if (!root)
        throw ReadError("Expected tag 'VTKFile': " + rPath);
    if (std::string(root.attribute("type").as_string()) != "vtkMultiBlockDataSet")
        throw ReadError("Expected type vtkMultiBlockDataSet: " + rPath);
    const pugi::xml_node mb = root.child("vtkMultiBlockDataSet");
    if (!mb)
        throw ReadError("Expected tag 'vtkMultiBlockDataSet': " + rPath);

    std::vector<pugi::xml_node> datasets;
    vtm_collect_datasets(mb, datasets);

    const fs::path base = fs::path(rPath).parent_path();
    std::vector<vtm_piece_ref> pieces;
    pieces.reserve(datasets.size());
    for (const pugi::xml_node& ds : datasets) {
        const std::string file = ds.attribute("file").as_string("");
        if (file.empty())
            throw ReadError("<DataSet> is missing its 'file' attribute: " + rPath);
        vtm_piece_ref ref;
        ref.mPath = base.empty() ? fs::path(file) : base / file;
        ref.mName = ds.attribute("name").as_string("");
        if (ref.mName.empty())
            ref.mName = "block_" + std::to_string(pieces.size());
        pieces.push_back(std::move(ref));
    }
    return pieces;
}

}  // namespace

void write_vtm_codec(const std::string& rPath, const Mesh& rMesh, bool binary, detail::VtkCodec codec) {
    if (binary && codec != detail::VtkCodec::None)
        detail::vtk_codec_require_write(codec);

    const fs::path index_path(rPath);
    const std::string stem = index_path.stem().string();
    const fs::path dir = index_path.parent_path().empty() ? fs::path(stem)
                                                            : index_path.parent_path() / stem;
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec)
        throw WriteError("Could not create directory for .vtm pieces: " + dir.string() + ": " +
                          ec.message());

    std::vector<std::string> piece_files;
    std::vector<std::string> piece_names;
    piece_files.reserve(rMesh.NumCellBlocks());
    piece_names.reserve(rMesh.NumCellBlocks());

    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        const Mesh piece = vtm_extract_piece(rMesh, i);
        CleanOptions copts;
        copts.weld = false;
        copts.remove_orphans = true;
        copts.drop_degenerate = false;
        copts.drop_duplicate_cells = false;
        const CleanResult cleaned = clean(piece, copts);

        const std::string piece_file_name = stem + "_" + std::to_string(i) + ".vtu";
        const fs::path piece_path = dir / piece_file_name;
        write_vtu_codec(piece_path.string(), cleaned.mMesh, binary, binary ? codec : detail::VtkCodec::None);

        piece_files.push_back((fs::path(stem) / piece_file_name).generic_string());
        piece_names.push_back("block_" + std::to_string(i));
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<vtkMultiBlockDataSet>\n";
    os << "<Block index=\"0\">\n";
    for (std::size_t i = 0; i < piece_files.size(); ++i)
        os << "<DataSet index=\"" << i << "\" name=\"" << piece_names[i] << "\" file=\""
           << piece_files[i] << "\"/>\n";
    os << "</Block>\n";
    os << "</vtkMultiBlockDataSet>\n";
    os << "</VTKFile>\n";
    if (!os)
        throw WriteError("Failed while writing: " + rPath);
}

void write_vtm(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_vtm_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

Mesh read_vtm(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    const std::vector<vtm_piece_ref> refs = vtm_parse_index(rPath, doc);
    if (refs.empty()) {
        Mesh empty;
        empty.AssignPoints(NDArray::Uninit(DType::Float64, {0, 3}));
        return empty;
    }

    std::vector<Mesh> pieces;
    pieces.reserve(refs.size());
    for (const vtm_piece_ref& ref : refs) {
        const std::string ext = ref.mPath.extension().string();
        if (ext == ".vtu")
            pieces.push_back(read_vtu(ref.mPath.string(), rOpts));
        else if (ext == ".vtp")
            pieces.push_back(read_vtp(ref.mPath.string(), rOpts));
        else
            throw ReadError("Unsupported .vtm piece '" + ref.mPath.string() +
                             "': only .vtu/.vtp pieces are read");
    }

    std::vector<const Mesh*> ptrs;
    ptrs.reserve(pieces.size());
    for (const Mesh& m : pieces)
        ptrs.push_back(&m);

    MergeOptions mopts;
    mopts.weld = false;
    mopts.source_tag = true;
    mopts.data_policy = MergeDataPolicy::Fill;
    MergeResult result = merge(ptrs, mopts);

    for (std::size_t i = 0; i < refs.size(); ++i)
        result.mMesh.AddRegion(Region(refs[i].mName, RegionKind::Cell, std::move(result.mCellMaps[i])));

    return std::move(result.mMesh);
}

MeshMetadata read_vtm_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    const std::vector<vtm_piece_ref> refs = vtm_parse_index(rPath, doc);

    MeshMetadata meta;
    meta.mFormat = "vtm";
    if (refs.empty())
        return meta;

    std::vector<CellBlockInfo> blocks;
    std::unordered_map<std::string, std::size_t> type_to_idx;
    std::vector<std::string> point_names, cell_names, field_names;
    auto merge_names = [](std::vector<std::string>& rInto, const std::vector<std::string>& rFrom) {
        for (const std::string& n : rFrom)
            if (std::find(rInto.begin(), rInto.end(), n) == rInto.end())
                rInto.push_back(n);
    };

    for (const vtm_piece_ref& ref : refs) {
        const std::string ext = ref.mPath.extension().string();
        const MeshMetadata pm = (ext == ".vtp") ? read_vtp_metadata(ref.mPath.string(), rOpts)
                                                 : read_vtu_metadata(ref.mPath.string(), rOpts);
        meta.mNumPoints += pm.mNumPoints;
        meta.mPointDim = std::max(meta.mPointDim, pm.mPointDim);
        meta.mFellBackToFullRead = meta.mFellBackToFullRead || pm.mFellBackToFullRead;
        merge_names(point_names, pm.mPointDataNames);
        merge_names(cell_names, pm.mCellDataNames);
        merge_names(field_names, pm.mFieldDataNames);

        for (const CellBlockInfo& cb : pm.mCellBlocks) {
            auto it = type_to_idx.find(cb.mType);
            if (it == type_to_idx.end()) {
                type_to_idx[cb.mType] = blocks.size();
                blocks.push_back(cb);
            } else {
                CellBlockInfo& dst = blocks[it->second];
                dst.mNumCells += cb.mNumCells;
                dst.mRagged = dst.mRagged || cb.mRagged;
                if (dst.mNodesPerCell != cb.mNodesPerCell)
                    dst.mNodesPerCell = 0;
            }
        }
    }

    meta.mCellBlocks = std::move(blocks);
    std::sort(point_names.begin(), point_names.end());
    std::sort(cell_names.begin(), cell_names.end());
    std::sort(field_names.begin(), field_names.end());
    meta.mPointDataNames = std::move(point_names);
    meta.mCellDataNames = std::move(cell_names);
    meta.mFieldDataNames = std::move(field_names);
    return meta;
}

}  // namespace meshioplusplus
