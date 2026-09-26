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
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/vtk_common.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/row_writer.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::read_double;
using detail::read_int;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

}  // namespace

void write_vtu(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    // The historical bool API, preserved exactly: zlib stays the only codec it
    // can select, so existing callers are unaffected.
    write_vtu_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_vtu_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t pt_isz = dtype_size(points.Dtype());

    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();

    const char* fmt = binary ? "binary" : "ascii";
    // UInt64 size headers only where an uncompressed array could pass 4 GiB
    // (compressed ones count 32 KiB blocks); everything else keeps its bytes.
    const std::size_t hsz =
        (binary && codec == detail::VtkCodec::None) ? detail::vtk_xml_header_bytes(rMesh) : 4;

    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0)
            os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\">\n";
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        os << detail::vtu_encode_binary(d, n, binary ? codec : detail::VtkCodec::None, hsz) << "\n";
    };

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    if (hsz == 8)
        os << " header_type=\"UInt64\"";
    os << ">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<UnstructuredGrid>\n";
    // Field data belongs to the dataset, not to a piece: VTK writes it on the grid,
    // before the <Piece>. Guarded, so a mesh without any writes the bytes it always did.
    if (rMesh.NumFieldData() != 0) {
        os << "<FieldData>\n";
        for (const auto& name : rMesh.FieldDataNames())
            detail::vtu_write_field_array(os, name, rMesh.FieldData(name), binary,
                                          binary ? codec : detail::VtkCodec::None, hsz);
        os << "</FieldData>\n";
    }
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << total_cells
       << "\">\n";

    // Points (3 components; pad 2D with zero z).
    os << "<Points>\n";
    da_header(vtu_type_str(points.Dtype()), "Points", 3);
    if (binary) {
        // Pre-sized buffer (zero-filled -> the padded z stays 0), indexed
        // byte writes -> parallel over points.
        std::vector<unsigned char> buf(num_points * 3 * pt_isz, 0);
        const auto* src = reinterpret_cast<const unsigned char*>(points.Data());
        parallel_for(num_points, [&](std::size_t r) {
            for (std::size_t c = 0; c < dim && c < 3; ++c)
                std::memcpy(buf.data() + (r * 3 + c) * pt_isz, src + (r * dim + c) * pt_isz,
                            pt_isz);
        });
        emit_bin(buf.data(), buf.size());
    } else {
        const detail::DoubleView pv(points);
        const detail::CNumber num;
        detail::write_row_chunks(
            os, num_points, [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                for (std::size_t r = First; r < Last; ++r)
                    for (std::size_t c = 0; c < 3; ++c) {
                        num.Append(rBuf, "%.11e", (c < dim) ? pv[r * dim + c] : 0.0);
                        rBuf += '\n';
                    }
            });
    }
    os << "</DataArray>\n</Points>\n";

    if (rMesh.NumCellBlocks() != 0) {
        // Build connectivity / offsets / types (Int64) into pre-sized arrays.
        // Per-block offsets are closed-form (conn_base + (r+1)*k), so rows are
        // independent and each block fills in parallel.
        const auto& tmap = meshio_to_vtk_type();
        std::size_t ncells = 0;
        bool any_polyhedron = false;
        for (const auto cb : rMesh.CellRange()) {
            ncells += cb.NumCells();
            if (cb.IsPolyhedron())
                any_polyhedron = true;
        }
        std::vector<std::int64_t> connectivity, offsets, types;
        offsets.reserve(ncells);
        types.reserve(ncells);
        // VTK_POLYHEDRON's face stream, and one END offset per cell into it.
        // A NON-polyhedral cell carries -1 -- that is exactly how VTU expresses
        // a mesh mixing polyhedra with other types, which is why meshio++
        // supports mixing where the Python reference historically did not (an
        // OpenFOAM mesh always mixes hexahedra, polyhedra and boundary faces).
        std::vector<std::int64_t> faces, face_offsets;
        if (any_polyhedron)
            face_offsets.reserve(ncells);

        for (const auto cb : rMesh.CellRange()) {
            const std::size_t nc = cb.NumCells();
            if (cb.IsPolyhedron()) {
                for (std::size_t r = 0; r < nc; ++r) {
                    // The cell-node connectivity VTK wants is the SORTED UNIQUE
                    // node set of the cell; the face structure lives entirely
                    // in `faces`. Sorted-unique matches the Python reference
                    // exactly -- a different order would give ParaView a
                    // different point set for the same cell.
                    std::vector<std::int64_t> uniq;
                    faces.push_back(static_cast<std::int64_t>(cb.NumFaces(r)));
                    for (std::size_t f = 0; f < cb.NumFaces(r); ++f) {
                        const auto face = cb.Face(r, f);
                        faces.push_back(static_cast<std::int64_t>(face.second));
                        for (std::size_t j = 0; j < face.second; ++j) {
                            faces.push_back(face.first[j]);
                            uniq.push_back(face.first[j]);
                        }
                    }
                    std::sort(uniq.begin(), uniq.end());
                    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                    connectivity.insert(connectivity.end(), uniq.begin(), uniq.end());
                    offsets.push_back(static_cast<std::int64_t>(connectivity.size()));
                    types.push_back(42);
                    face_offsets.push_back(static_cast<std::int64_t>(faces.size()));
                }
                continue;
            }
            if (cb.IsRagged()) {
                // A jagged polygon block: VTK_POLYGON with a per-cell size.
                for (std::size_t r = 0; r < nc; ++r) {
                    for (std::size_t j = 0; j < cb.RowSize(r); ++j)
                        connectivity.push_back(cb.Row(r)[j]);
                    offsets.push_back(static_cast<std::int64_t>(connectivity.size()));
                    types.push_back(7);  // VTK_POLYGON
                    if (any_polyhedron)
                        face_offsets.push_back(-1);
                }
                continue;
            }
            const NDArray& conn = cb.Conn();
            const std::size_t k = cols(conn);
            std::vector<int> order = meshio_to_vtk_order(cb.Type());
            auto it = tmap.find(cb.Type());
            if (it == tmap.end())
                throw WriteError("Unknown cell type for VTU: " + cb.Type());
            const std::int64_t vtk_type = it->second;
            // Rows are independent: each one's slots and offset are closed-form.
            const detail::Int64View cv(conn);
            const std::size_t cbase = connectivity.size();
            const std::size_t rbase = offsets.size();
            connectivity.resize(cbase + nc * k);
            offsets.resize(rbase + nc);
            types.resize(rbase + nc, vtk_type);
            if (any_polyhedron)
                face_offsets.resize(face_offsets.size() + nc, -1);
            parallel_for_bw(nc, [&](std::size_t r) {
                for (std::size_t j = 0; j < k; ++j) {
                    const std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    connectivity[cbase + r * k + j] = cv[r * k + col];
                }
                offsets[rbase + r] = static_cast<std::int64_t>(cbase + (r + 1) * k);
            });
        }

        auto emit_i64 = [&](const char* name, const std::vector<std::int64_t>& v) {
            da_header("Int64", name, 0);
            if (binary) {
                emit_bin(reinterpret_cast<const unsigned char*>(v.data()),
                         v.size() * sizeof(std::int64_t));
            } else {
                detail::write_row_chunks(
                    os, v.size(), [&](std::size_t First, std::size_t Last, std::string& rBuf) {
                        for (std::size_t i = First; i < Last; ++i) {
                            detail::append_int(rBuf, v[i]);
                            rBuf += '\n';
                        }
                    });
            }
            os << "</DataArray>\n";
        };

        os << "<Cells>\n";
        emit_i64("connectivity", connectivity);
        emit_i64("offsets", offsets);
        emit_i64("types", types);
        if (any_polyhedron) {
            emit_i64("faces", faces);
            emit_i64("faceoffsets", face_offsets);
        }
        os << "</Cells>\n";
    }

    if (rMesh.NumPointData() != 0) {
        os << "<PointData>\n";
        for (const auto& name : rMesh.PointDataNames()) {
            NDArray scratch;
            const NDArray& d = detail::vtu_disk_array(name, rMesh.PointData(name), scratch);
            int ncomp = (d.Shape().size() == 2) ? static_cast<int>(cols(d)) : 0;
            da_header(vtu_type_str(d.Dtype()), name, ncomp);
            if (binary)
                emit_bin(reinterpret_cast<const unsigned char*>(d.Data()), d.Nbytes());
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
            NDArray scratch;
            const NDArray& first = rMesh.CellData(name, 0);
            int ncomp = (first.Shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            da_header(vtu_type_str(detail::vtu_disk_dtype(name, first.Dtype())), name, ncomp);
            if (binary && nblocks == 1) {
                // One block: its bytes go out as they are, with no copy.
                const NDArray& blk = detail::vtu_disk_array(name, first, scratch);
                emit_bin(reinterpret_cast<const unsigned char*>(blk.Data()), blk.Nbytes());
            } else if (binary) {
                std::size_t total = 0;
                const std::size_t isz = dtype_size(detail::vtu_disk_dtype(name, first.Dtype()));
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    total += rMesh.CellData(name, bi).Size() * isz;
                std::vector<unsigned char> buf;
                buf.reserve(total);
                for (std::size_t bi = 0; bi < nblocks; ++bi) {
                    const NDArray& blk =
                        detail::vtu_disk_array(name, rMesh.CellData(name, bi), scratch);
                    const unsigned char* p = reinterpret_cast<const unsigned char*>(blk.Data());
                    buf.insert(buf.end(), p, p + blk.Nbytes());
                }
                emit_bin(buf.data(), buf.size());
            } else {
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    vtu_ascii_ndarray(
                        os, detail::vtu_disk_array(name, rMesh.CellData(name, bi), scratch));
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace meshioplusplus
