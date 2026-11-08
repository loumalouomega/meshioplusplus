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
 * @file detail/vtk_xml_write.hpp
 * @brief The pieces the VTK XML writers (`.vti`, `.vtr`, `.vts`, `.vtp`,
 * `.vtu`) share: the `<VTKFile>` prolog, the DataArray header and body
 * emitters, and a lattice's PointData/CellData sections.
 *
 * A **core-private** header: no installed header names it, and it adds
 * nothing to the API or the ABI.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "region_field_data.hpp"
#include "row_writer.hpp"
#include "vtu_decode.hpp"

namespace meshioplusplus {
namespace detail {

/// The binary size-header width: UInt64 size headers only where an
/// uncompressed array could pass 4 GiB (compressed ones count 32 KiB
/// blocks); everything else keeps its bytes.
inline std::size_t vtk_xml_size_header_bytes(const Mesh& rMesh, bool Binary, VtkCodec Codec) {
    return (Binary && Codec == VtkCodec::None) ? vtk_xml_header_bytes(rMesh) : 4;
}

/// The XML declaration, the `<VTKFile>` start tag and the provenance comment.
inline void vtk_xml_prolog(std::ostream& rOs, const char* pType, const char* pVersion, bool Binary,
                           VtkCodec Codec, std::size_t HeaderBytes) {
    rOs << "<?xml version=\"1.0\"?>\n";
    rOs << "<VTKFile type=\"" << pType << "\" version=\"" << pVersion
        << "\" byte_order=\"LittleEndian\"";
    if (Binary && Codec != VtkCodec::None)
        rOs << " compressor=\"" << vtk_codec_compressor(Codec) << "\"";
    if (HeaderBytes == 8)
        rOs << " header_type=\"UInt64\"";
    rOs << ">\n";
    rOs << provenance_render_xml_comment(SlotTier::Block) << "\n";
}

/// One Int64 DataArray through a writer's own header and binary emitters.
template <class Header, class Bin>
void vtk_xml_emit_i64(std::ostream& rOs, bool Binary, Header&& DaHeader, Bin&& EmitBin,
                      const char* pName, const std::vector<std::int64_t>& rV) {
    DaHeader("Int64", pName, 0);
    if (Binary) {
        EmitBin(reinterpret_cast<const unsigned char*>(rV.data()),
                rV.size() * sizeof(std::int64_t));
    } else {
        write_row_chunks(rOs, rV.size(),
                         [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                             for (std::size_t i = First; i < Last; ++i) {
                                 append_int(rBuf, rV[i]);
                                 rBuf += '\n';
                             }
                         });
    }
    rOs << "</DataArray>\n";
}

/// The inline (ascii or binary) DataArray emitters of a VTK XML writer.
struct VtkXmlEmitter {
    std::ostream& mOs;
    bool mBinary;
    VtkCodec mCodec;
    std::size_t mHeaderBytes;

    /// A `<DataArray>` start tag; @p Ncomp 0 leaves NumberOfComponents out.
    void Header(const char* pType, const std::string& rName, int Ncomp) const {
        mOs << "<DataArray type=\"" << pType << "\" Name=\"" << rName << "\"";
        if (Ncomp > 0)
            mOs << " NumberOfComponents=\"" << Ncomp << "\"";
        mOs << " format=\"" << (mBinary ? "binary" : "ascii") << "\">\n";
    }

    /// A binary DataArray body and its line end.
    void Bin(const unsigned char* pData, std::size_t N) const {
        vtu_write_binary_bytes(mOs, pData, N, mBinary ? mCodec : VtkCodec::None, mHeaderBytes);
        mOs << "\n";
    }

    /// One Int64 DataArray.
    void Int64s(const char* pName, const std::vector<std::int64_t>& rV) const {
        vtk_xml_emit_i64(
            mOs, mBinary,
            [this](const char* pType, const std::string& rName, int Ncomp) {
                Header(pType, rName, Ncomp);
            },
            [this](const unsigned char* pData, std::size_t N) { Bin(pData, N); }, pName, rV);
    }
};

/// A lattice's `WholeExtent`/`Extent` attribute: "0 nx 0 ny 0 nz".
inline std::string vtk_xml_lattice_extent(const LatticeSpec& rSpec) {
    auto ext = make_classic_ostringstream();
    ext << "0 " << rSpec.mDims[0] << " 0 " << rSpec.mDims[1] << " 0 " << rSpec.mDims[2];
    return ext.str();
}

/// The opening of a lattice writer whose grid element carries only the extent
/// (`.vtr`, `.vts`): the prolog, the `<pGrid WholeExtent="...">` start tag,
/// the field data section and the `<Piece>` start tag. @p pFormat names the
/// format in the field data warnings. Returns the size-header width.
inline std::size_t vtk_xml_lattice_preamble(std::ostream& rOs, const Mesh& rMesh,
                                            const LatticeSpec& rSpec, const char* pGrid,
                                            const char* pFormat, bool Binary, VtkCodec Codec) {
    const std::size_t hsz = vtk_xml_size_header_bytes(rMesh, Binary, Codec);
    const std::string ext = vtk_xml_lattice_extent(rSpec);

    vtk_xml_prolog(rOs, pGrid, "0.1", Binary, Codec, hsz);
    rOs << "<" << pGrid << " WholeExtent=\"" << ext << "\">\n";
    write_field_data_section(rOs, rMesh, pFormat, Binary, Codec, hsz);
    rOs << "<Piece Extent=\"" << ext << "\">\n";
    return hsz;
}

/// The codec of the `(binary, zlib)` writer entry points: zlib or none.
inline VtkCodec vtk_codec_from_zlib(bool Zlib) {
    return Zlib ? VtkCodec::Zlib : VtkCodec::None;
}

/**
 * @brief The start of every lattice writer (`.vti`, `.vtr`, `.vts`): the
 * mesh's lattice into @p rSpec, the codec checked, the file opened.
 * @throws WriteError @p pNotLattice when the mesh is not a lattice, and when
 * the file cannot be opened.
 */
inline auto vtk_xml_open_lattice(const std::string& rPath, const Mesh& rMesh, bool Binary,
                                 VtkCodec Codec, const char* pNotLattice,
                                 LatticeSpec& rSpec) -> std::ofstream {
    if (!lattice_from_mesh(rMesh, rSpec))
        throw WriteError(pNotLattice);
    if (Binary && Codec != VtkCodec::None)
        vtk_codec_require_write(Codec);
    auto os = make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    return os;
}

/// A one-block lattice's `<PointData>` and `<CellData>` sections (each only
/// when the mesh has such data).
inline void vtk_xml_lattice_data(const VtkXmlEmitter& rEmit, const Mesh& rMesh) {
    std::ostream& os = rEmit.mOs;
    if (rMesh.NumPointData() != 0) {
        os << "<PointData>\n";
        for (const auto& name : rMesh.PointDataNames()) {
            const NDArray& d = rMesh.PointData(name);
            const int ncomp = (d.Shape().size() == 2) ? static_cast<int>(cols(d)) : 0;
            rEmit.Header(vtu_type_str(d.Dtype()), name, ncomp);
            if (rEmit.mBinary)
                rEmit.Bin(reinterpret_cast<const unsigned char*>(d.Data()), d.Nbytes());
            else
                vtu_ascii_ndarray(os, d);
            os << "</DataArray>\n";
        }
        os << "</PointData>\n";
    }

    if (rMesh.NumCellData() != 0) {
        os << "<CellData>\n";
        for (const auto& name : rMesh.CellDataNames()) {
            const std::size_t nblocks = rMesh.CellDataNumBlocks(name);
            if (nblocks == 0)
                continue;
            // A lattice has exactly one block, so there is nothing to
            // concatenate -- but iterate anyway rather than assume, since an
            // array that does not cover the block is a caller error worth not
            // writing silently truncated.
            const NDArray& first = rMesh.CellData(name, 0);
            const int ncomp = (first.Shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            rEmit.Header(vtu_type_str(first.Dtype()), name, ncomp);
            if (rEmit.mBinary) {
                std::vector<unsigned char> buf;
                for (std::size_t bi = 0; bi < nblocks; ++bi) {
                    const NDArray& blk = rMesh.CellData(name, bi);
                    const auto* p = reinterpret_cast<const unsigned char*>(blk.Data());
                    buf.insert(buf.end(), p, p + blk.Nbytes());
                }
                rEmit.Bin(buf.data(), buf.size());
            } else {
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    vtu_ascii_ndarray(os, rMesh.CellData(name, bi));
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }
}

/**
 * @brief A whole `.vtr` or `.vts` file, whose grid element carries only the
 * extent: open (`vtk_xml_open_lattice`), the preamble, the geometry section
 * `rGeometry(emit, spec)` writes (`<Coordinates>`, `<Points>`), the data
 * sections and the closing tags.
 */
template <class TGeometry>
void vtk_xml_write_lattice(const std::string& rPath, const Mesh& rMesh, bool Binary, VtkCodec Codec,
                           const char* pGrid, const char* pFormat, const char* pNotLattice,
                           TGeometry&& rGeometry) {
    LatticeSpec spec;
    auto os = vtk_xml_open_lattice(rPath, rMesh, Binary, Codec, pNotLattice, spec);
    const std::size_t hsz =
        vtk_xml_lattice_preamble(os, rMesh, spec, pGrid, pFormat, Binary, Codec);
    const VtkXmlEmitter emit{os, Binary, Codec, hsz};
    rGeometry(emit, spec);
    vtk_xml_lattice_data(emit, rMesh);
    os << "</Piece>\n</" << pGrid << ">\n</VTKFile>\n";
}

}  // namespace detail
}  // namespace meshioplusplus
