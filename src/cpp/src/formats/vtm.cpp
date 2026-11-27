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
//                   SPDX-License-Identifier: MIT
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
#include "pindex_common.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/ragged_csr.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "../detail/region_tree.hpp"
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
    if (cb.IsRagged()) {
        detail::append_ragged_copy(cb, piece);
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

// Every `<DataSet>` under `rNode`, in document order, tracking ancestor
// `<Block name="...">` hierarchy prefixes and empty named `<Block>` groups.
bool vtm_collect_datasets(const pugi::xml_node& rNode, const std::string& rPrefix,
                          std::vector<std::pair<std::string, pugi::xml_node>>& rOut,
                          std::vector<std::string>& rEmptyGroups) {
    bool found_any = false;
    for (pugi::xml_node child : rNode.children()) {
        const std::string tag(child.name());
        if (tag == "DataSet") {
            rOut.emplace_back(rPrefix, child);
            found_any = true;
        } else if (tag == "Block") {
            const std::string bname = child.attribute("name").as_string("");
            const std::string next_prefix =
                bname.empty() ? rPrefix : (rPrefix.empty() ? bname : rPrefix + "/" + bname);
            const bool sub = vtm_collect_datasets(child, next_prefix, rOut, rEmptyGroups);
            if (!sub && !bname.empty())
                rEmptyGroups.push_back(next_prefix);
            found_any = found_any || sub;
        } else {
            if (vtm_collect_datasets(child, rPrefix, rOut, rEmptyGroups))
                found_any = true;
        }
    }
    return found_any;
}

// One piece, resolved from its `<DataSet>` element: the path to read (relative
// to the index file's own directory) and the name to give it.
struct vtm_piece_ref {
    fs::path mPath;
    std::string mName;
};

std::vector<vtm_piece_ref> vtm_parse_index(const std::string& rPath, pugi::xml_document& rDoc,
                                           std::vector<std::string>* pEmptyGroups = nullptr) {
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

    std::vector<std::pair<std::string, pugi::xml_node>> datasets;
    std::vector<std::string> empty_groups;
    vtm_collect_datasets(mb, "", datasets, empty_groups);
    if (pEmptyGroups)
        *pEmptyGroups = std::move(empty_groups);

    const fs::path base = fs::path(rPath).parent_path();
    std::vector<vtm_piece_ref> pieces;
    pieces.reserve(datasets.size());
    for (const auto& [prefix, ds] : datasets) {
        const std::string file = ds.attribute("file").as_string("");
        if (file.empty())
            throw ReadError("<DataSet> is missing its 'file' attribute: " + rPath);
        vtm_piece_ref ref;
        ref.mPath = base.empty() ? fs::path(file) : base / file;
        std::string leaf = ds.attribute("name").as_string("");
        if (leaf.empty())
            leaf = "block_" + std::to_string(pieces.size());
        ref.mName = prefix.empty() ? leaf : (prefix + "/" + leaf);
        pieces.push_back(std::move(ref));
    }
    return pieces;
}

// A piece's name: the cell region that is exactly that block's cells (so a mesh
// read from a .vtm writes its pieces back under their names), else block_<i>.
// A name used by an earlier piece, or none, falls back to block_<i>.
std::vector<std::string> vtm_piece_names(const Mesh& rMesh) {
    std::vector<std::int64_t> starts(rMesh.NumCellBlocks() + 1, 0);
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i)
        starts[i + 1] = starts[i] + static_cast<std::int64_t>(rMesh.Cells(i).NumCells());
    std::vector<std::string> names;
    for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i) {
        std::string name = "block_" + std::to_string(i);
        const std::size_t count = static_cast<std::size_t>(starts[i + 1] - starts[i]);
        for (std::size_t r = 0; r < rMesh.NumRegions() && count > 0; ++r) {
            const Region& reg = rMesh.Region(r);
            if (reg.mKind != RegionKind::Cell || reg.mName.empty() ||
                reg.NumEntries() != count || std::find(names.begin(), names.end(), reg.mName) != names.end())
                continue;
            const std::int64_t* e = reg.Entries();
            bool exact = true;
            for (std::size_t k = 0; k < count && exact; ++k)
                exact = e[k] == starts[i] + static_cast<std::int64_t>(k);
            if (exact) {
                name = reg.mName;
                break;
            }
        }
        names.push_back(std::move(name));
    }
    return names;
}

std::string vtm_xml_attr(const std::string& rText) {
    std::string out;
    for (const char c : rText) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c;
        }
    }
    return out;
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
    const std::vector<std::string> piece_names = vtm_piece_names(rMesh);
    piece_files.reserve(rMesh.NumCellBlocks());

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
    }

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<vtkMultiBlockDataSet>\n";

    const std::vector<std::string> empty_tree_nodes = detail::get_region_tree_nodes(rMesh);
    bool has_hierarchy = !empty_tree_nodes.empty();
    for (const auto& name : piece_names)
        if (name.find('/') != std::string::npos)
            has_hierarchy = true;

    if (!has_hierarchy) {
        os << "<Block index=\"0\">\n";
        for (std::size_t i = 0; i < piece_files.size(); ++i)
            os << "<DataSet index=\"" << i << "\" name=\"" << vtm_xml_attr(piece_names[i])
               << "\" file=\"" << piece_files[i] << "\"/>\n";
        os << "</Block>\n";
    } else {
        struct VtmTreeNode {
            std::string mName;
            std::vector<std::pair<std::string, std::size_t>> mDatasets;  // (leaf_name, piece_idx)
            std::vector<std::string> mChildOrder;
            std::unordered_map<std::string, VtmTreeNode> mChildren;
        };
        VtmTreeNode root_node;
        auto ensure_group = [&](const std::vector<std::string>& rSegs,
                                std::size_t Depth) -> VtmTreeNode* {
            VtmTreeNode* cur = &root_node;
            for (std::size_t d = 0; d < Depth; ++d) {
                const std::string& seg = rSegs[d];
                auto it = cur->mChildren.find(seg);
                if (it == cur->mChildren.end()) {
                    cur->mChildOrder.push_back(seg);
                    VtmTreeNode child;
                    child.mName = seg;
                    it = cur->mChildren.emplace(seg, std::move(child)).first;
                }
                cur = &it->second;
            }
            return cur;
        };
        for (std::size_t i = 0; i < piece_names.size(); ++i) {
            const auto segs = detail::split_region_path(piece_names[i]);
            if (segs.empty()) {
                root_node.mDatasets.emplace_back("block_" + std::to_string(i), i);
            } else {
                VtmTreeNode* parent = ensure_group(segs, segs.size() - 1);
                parent->mDatasets.emplace_back(segs.back(), i);
            }
        }
        for (const std::string& p : empty_tree_nodes) {
            const auto segs = detail::split_region_path(p);
            if (!segs.empty())
                (void)ensure_group(segs, segs.size());
        }
        auto emit_node = [&](auto&& self, const VtmTreeNode& rNode) -> void {
            std::size_t idx = 0;
            for (const auto& [leaf, pi] : rNode.mDatasets) {
                os << "<DataSet index=\"" << idx++ << "\" name=\"" << vtm_xml_attr(leaf)
                   << "\" file=\"" << piece_files[pi] << "\"/>\n";
            }
            for (const std::string& cname : rNode.mChildOrder) {
                const VtmTreeNode& child = rNode.mChildren.at(cname);
                os << "<Block index=\"" << idx++ << "\" name=\"" << vtm_xml_attr(child.mName)
                   << "\">\n";
                self(self, child);
                os << "</Block>\n";
            }
        };
        emit_node(emit_node, root_node);
    }

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
    std::vector<std::string> empty_groups;
    const std::vector<vtm_piece_ref> refs = vtm_parse_index(rPath, doc, &empty_groups);
    if (refs.empty()) {
        Mesh empty;
        empty.AssignPoints(NDArray::Uninit(DType::Float64, {0, 3}));
        if (!empty_groups.empty())
            detail::set_region_tree_nodes(empty, empty_groups);
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
    MergeResult result = pidx::merge_keeping_field_data(ptrs, mopts);

    for (std::size_t i = 0; i < refs.size(); ++i)
        result.mMesh.AddRegion(Region(refs[i].mName, RegionKind::Cell, std::move(result.mCellMaps[i])));
    if (!empty_groups.empty())
        detail::set_region_tree_nodes(result.mMesh, empty_groups);

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
