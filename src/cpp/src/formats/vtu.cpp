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
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtu.hpp"
#include "../detail/region_field_data.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/vtk_common.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "../detail/row_writer.hpp"
#include "../detail/typed_view.hpp"
#include "../detail/vtu_decode.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::read_double;
using detail::read_int;
using detail::vtu_ascii_ndarray;
using detail::vtu_type_str;

// Bytes [Offset, Offset + N) of an array of `Isz`-byte items into pDst, for
// the streaming encoders (detail::VtuFill): `rItems(First, Last, pOut)` writes
// the items [First, Last) to pOut. The encoders ask for whole 32 KiB blocks
// and every item is 1, 2, 4 or 8 bytes, so a window never cuts an item.
template <class F>
void vtu_fill_items(std::size_t Offset, std::size_t N, unsigned char* pDst, std::size_t Isz,
                    const F& rItems) {
    if (Offset % Isz != 0 || N % Isz != 0)
        throw WriteError("VTU: a streamed window cuts an item of " + std::to_string(Isz) +
                         " bytes");
    rItems(Offset / Isz, (Offset + N) / Isz, pDst);
}

// One fixed-size cell block of a VTU piece: where its cells and connectivity
// start in the piece's arrays, and how a row maps to VTK's node order.
struct VtuBlockRun {
    std::size_t mCellBase = 0;
    std::size_t mConnBase = 0;
    std::size_t mNumCells = 0;
    std::size_t mK = 0;
    std::int64_t mType = 0;
    std::vector<int> mOrder;
    std::unique_ptr<detail::Int64View> mpConn;
};

// The int64 items [First, Last) of a piece-wide array whose items belong to
// blocks by the range [Base(b), Base(b) + Count(b)), as `rValue(b, i)` (i local
// to the block), written to pOut -- each block's share in parallel.
template <class B, class C, class V>
void vtu_fill_block_items(const std::vector<VtuBlockRun>& rRuns, std::size_t First,
                          std::size_t Last, unsigned char* pOut, const B& rBase, const C& rCount,
                          const V& rValue) {
    for (const VtuBlockRun& run : rRuns) {
        const std::size_t base = rBase(run);
        const std::size_t lo = std::max(First, base);
        const std::size_t hi = std::min(Last, base + rCount(run));
        if (lo >= hi)
            continue;
        parallel_for_bw(hi - lo, [&](std::size_t t) {
            const std::int64_t v = rValue(run, lo + t - base);
            std::memcpy(pOut + (lo + t - First) * sizeof(std::int64_t), &v, sizeof(v));
        });
    }
}

// The writer behind write_vtu_codec and write_vtu_appended. `Appended` writes
// every array as `format="appended"`: its size header and payload (or
// block header and compressed blocks) go, unencoded, into one raw
// `<AppendedData>` section at the end of the file, each array's `offset`
// counted from the byte after its leading underscore.
void vtu_write_impl(const std::string& rPath, const Mesh& rMesh, bool binary,
                    detail::VtkCodec codec, bool Appended) {
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

    const char* fmt = Appended ? "appended" : (binary ? "binary" : "ascii");
    std::vector<unsigned char> appended;  // the raw <AppendedData> payload
    // UInt64 size headers only where an uncompressed array could pass 4 GiB
    // (compressed ones count 32 KiB blocks); everything else keeps its bytes.
    const std::size_t hsz =
        (binary && codec == detail::VtkCodec::None) ? detail::vtk_xml_header_bytes(rMesh) : 4;

    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0)
            os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\"";
        if (Appended)
            os << " offset=\"" << appended.size() << "\"";
        os << ">\n";
    };
    // Every binary array goes out through the streaming encoders, which ask
    // `rFill` for its bytes a fixed batch at a time: an array generated on
    // the fly (the padded points, the cell arrays) is never held whole, nor
    // are its compressed blocks' copies or its base64 text (roadmap §6.1).
    auto emit_stream = [&](std::size_t n, const detail::VtuFill& rFill) {
        if (Appended) {
            detail::vtu_encode_raw_stream(n, codec, hsz, rFill, appended);
        } else {
            detail::vtu_write_binary_stream(os, n, binary ? codec : detail::VtkCodec::None, hsz,
                                            rFill);
            os << "\n";
        }
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        emit_stream(n, [d](std::size_t Off, std::size_t N, unsigned char* pDst) {
            std::memcpy(pDst, d + Off, N);
        });
    };

    // VTK XML 2.2 (VTK 9.1) is the first version whose higher-order hexahedra
    // list edges 10 and 11 in the order meshio++ holds; an older one tells VTK to
    // swap them (detail/vtk_hex_order.hpp), so a mesh with such cells says 2.2.
    bool hex_order = false;
    for (const auto cb : rMesh.CellRange())
        hex_order = hex_order || cb.Type() == "VTK_LAGRANGE_HEXAHEDRON" ||
                    cb.Type() == "VTK_BEZIER_HEXAHEDRON";
    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"" << (hex_order ? "2.2" : "0.1")
       << "\" byte_order=\"LittleEndian\"";
    if (binary && codec != detail::VtkCodec::None)
        os << " compressor=\"" << detail::vtk_codec_compressor(codec) << "\"";
    if (hsz == 8)
        os << " header_type=\"UInt64\"";
    os << ">\n";
    os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
    os << "<UnstructuredGrid>\n";
    // Field data belongs to the dataset, not to a piece: VTK writes it on the grid,
    // before the <Piece>. Guarded, so a mesh without any writes the bytes it always did.
    // Named regions ride here too (detail/region_field_data.hpp); a VTU's file
    // cell order is block-major, so no cell translation is needed.
    const std::vector<std::pair<std::string, NDArray>> region_arrays =
        detail::regions_to_field_arrays(rMesh, nullptr);
    std::vector<std::pair<std::string, const NDArray*>> field_arrays;
    for (const auto& name : rMesh.FieldDataNames()) {
        if (detail::is_region_field_name(name)) {
            log::warn("vtu: field_data '{}' uses the region naming convention; not written", name);
            continue;
        }
        const auto& array = rMesh.FieldData(name);
        const auto& shape = array.Shape();
        if (shape.size() > 1 && std::find(shape.begin() + 1, shape.end(), 0) != shape.end()) {
            log::warn("vtu: field_data '{}' has zero components; not written", name);
            continue;
        }
        field_arrays.emplace_back(name, &array);
    }
    for (const auto& [name, arr] : region_arrays)
        field_arrays.emplace_back(name, &arr);
    if (!field_arrays.empty()) {
        os << "<FieldData>\n";
        for (const auto& [name, p_arr] : field_arrays) {
            const NDArray& arr = *p_arr;
            if (!Appended) {
                detail::vtu_write_field_array(os, name, arr, binary,
                                              binary ? codec : detail::VtkCodec::None, hsz);
                continue;
            }
            // vtu_write_field_array's element, with the payload appended.
            const std::vector<std::size_t>& shape = arr.Shape();
            os << "<DataArray type=\"" << vtu_type_str(arr.Dtype()) << "\" Name=\"" << name
               << "\" NumberOfTuples=\"" << (shape.empty() ? 1 : shape[0]) << "\"";
            if (shape.size() >= 2) {
                std::size_t components = 1;
                for (std::size_t d = 1; d < shape.size(); ++d)
                    components *= shape[d];
                os << " NumberOfComponents=\"" << components << "\"";
            }
            os << " format=\"appended\" offset=\"" << appended.size() << "\">\n";
            emit_bin(reinterpret_cast<const unsigned char*>(arr.Data()), arr.Nbytes());
            os << "</DataArray>\n";
        }
        os << "</FieldData>\n";
    }
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << total_cells
       << "\">\n";

    // Points (3 components; pad 2D with zero z).
    os << "<Points>\n";
    da_header(vtu_type_str(points.Dtype()), "Points", 3);
    if (binary) {
        // Three components per point, generated a window at a time: a 3-D
        // array streams as it is, a 2-D one with a zero z per point.
        const auto* src = reinterpret_cast<const unsigned char*>(points.Data());
        if (dim == 3) {
            emit_bin(src, num_points * 3 * pt_isz);
        } else {
            emit_stream(
                num_points * 3 * pt_isz, [&](std::size_t Off, std::size_t N, unsigned char* pDst) {
                    vtu_fill_items(Off, N, pDst, pt_isz,
                                   [&](std::size_t First, std::size_t Last, unsigned char* pOut) {
                                       parallel_for(Last - First, [&](std::size_t t) {
                                           const std::size_t i = First + t;
                                           const std::size_t r = i / 3, c = i % 3;
                                           unsigned char* o = pOut + t * pt_isz;
                                           if (c < dim)
                                               std::memcpy(o, src + (r * dim + c) * pt_isz, pt_isz);
                                           else
                                               std::memset(o, 0, pt_isz);
                                       });
                                   });
                });
        }
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
        bool any_polyhedron = false, any_ragged = false;
        for (const auto cb : rMesh.CellRange()) {
            ncells += cb.NumCells();
            if (cb.IsPolyhedron())
                any_polyhedron = true;
            else if (cb.IsRagged())
                any_ragged = true;
        }
        if (binary && !any_polyhedron && !any_ragged) {
            // Fixed-size blocks only: every item of connectivity, offsets and
            // types is closed-form in its block, so the three arrays are
            // generated a window at a time by the streaming encoder rather
            // than built whole (they are 3 x 8 bytes per cell plus 8 per node).
            std::vector<VtuBlockRun> runs;
            std::size_t cell_base = 0, conn_base = 0;
            for (const auto cb : rMesh.CellRange()) {
                auto it = tmap.find(cb.Type());
                if (it == tmap.end())
                    throw WriteError("Unknown cell type for VTU: " + cb.Type());
                VtuBlockRun run;
                run.mCellBase = cell_base;
                run.mConnBase = conn_base;
                run.mNumCells = cb.NumCells();
                run.mK = cols(cb.Conn());
                run.mType = it->second;
                run.mOrder = meshio_to_vtk_order(cb.Type());
                run.mpConn = std::make_unique<detail::Int64View>(cb.Conn());
                cell_base += run.mNumCells;
                conn_base += run.mNumCells * run.mK;
                runs.push_back(std::move(run));
            }
            const auto cell_of = [](const VtuBlockRun& rRun) { return rRun.mCellBase; };
            const auto cells_in = [](const VtuBlockRun& rRun) { return rRun.mNumCells; };
            const auto stream_i64 = [&](const char* name, std::size_t count, const auto& rItems) {
                da_header("Int64", name, 0);
                emit_stream(count * sizeof(std::int64_t),
                            [&](std::size_t Off, std::size_t N, unsigned char* pDst) {
                                vtu_fill_items(Off, N, pDst, sizeof(std::int64_t), rItems);
                            });
                os << "</DataArray>\n";
            };
            os << "<Cells>\n";
            stream_i64(
                "connectivity", conn_base,
                [&](std::size_t First, std::size_t Last, unsigned char* pOut) {
                    vtu_fill_block_items(
                        runs, First, Last, pOut,
                        [](const VtuBlockRun& rRun) { return rRun.mConnBase; },
                        [](const VtuBlockRun& rRun) { return rRun.mNumCells * rRun.mK; },
                        [](const VtuBlockRun& rRun, std::size_t i) {
                            const std::size_t r = i / rRun.mK, j = i % rRun.mK;
                            const std::size_t col =
                                rRun.mOrder.empty() ? j : static_cast<std::size_t>(rRun.mOrder[j]);
                            return (*rRun.mpConn)[r * rRun.mK + col];
                        });
                });
            stream_i64("offsets", cell_base,
                       [&](std::size_t First, std::size_t Last, unsigned char* pOut) {
                           vtu_fill_block_items(runs, First, Last, pOut, cell_of, cells_in,
                                                [](const VtuBlockRun& rRun, std::size_t r) {
                                                    return static_cast<std::int64_t>(
                                                        rRun.mConnBase + (r + 1) * rRun.mK);
                                                });
                       });
            stream_i64("types", cell_base,
                       [&](std::size_t First, std::size_t Last, unsigned char* pOut) {
                           vtu_fill_block_items(
                               runs, First, Last, pOut, cell_of, cells_in,
                               [](const VtuBlockRun& rRun, std::size_t) { return rRun.mType; });
                       });
            os << "</Cells>\n";
        } else {
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
                        const std::size_t col =
                            order.empty() ? j : static_cast<std::size_t>(order[j]);
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
    }

    // A `RationalWeights` point array and a `HigherOrderDegrees` cell array are
    // marked as such (since v16.40.0), which is how VTK finds the weights and
    // degrees of Bezier cells; a mesh without them writes the same bytes.
    if (rMesh.NumPointData() != 0) {
        os << (rMesh.HasPointData("RationalWeights")
                   ? "<PointData RationalWeights=\"RationalWeights\">\n"
                   : "<PointData>\n");
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
        os << (rMesh.HasCellData("HigherOrderDegrees")
                   ? "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n"
                   : "<CellData>\n");
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
                // The blocks' bytes in turn, streamed from each block in
                // place rather than concatenated into one buffer first.
                std::vector<NDArray> scratches(nblocks);
                std::vector<const NDArray*> blks(nblocks);
                std::vector<std::size_t> starts(nblocks + 1, 0);
                for (std::size_t bi = 0; bi < nblocks; ++bi) {
                    blks[bi] =
                        &detail::vtu_disk_array(name, rMesh.CellData(name, bi), scratches[bi]);
                    starts[bi + 1] = starts[bi] + blks[bi]->Nbytes();
                }
                emit_stream(starts[nblocks], [&](std::size_t Off, std::size_t N,
                                                 unsigned char* pDst) {
                    for (std::size_t bi = 0; bi < nblocks; ++bi) {
                        const std::size_t lo = std::max(Off, starts[bi]);
                        const std::size_t hi = std::min(Off + N, starts[bi + 1]);
                        if (lo < hi)
                            std::memcpy(pDst + (lo - Off),
                                        reinterpret_cast<const unsigned char*>(blks[bi]->Data()) +
                                            (lo - starts[bi]),
                                        hi - lo);
                    }
                });
            } else {
                for (std::size_t bi = 0; bi < nblocks; ++bi)
                    vtu_ascii_ndarray(
                        os, detail::vtu_disk_array(name, rMesh.CellData(name, bi), scratch));
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</UnstructuredGrid>\n";
    if (Appended) {
        os << "<AppendedData encoding=\"raw\">\n_";
        os.write(reinterpret_cast<const char*>(appended.data()),
                 static_cast<std::streamsize>(appended.size()));
        os << "\n</AppendedData>\n";
    }
    os << "</VTKFile>\n";
}

}  // namespace

void write_vtu(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    // The historical bool API, preserved exactly: zlib stays the only codec it
    // can select, so existing callers are unaffected.
    write_vtu_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_vtu_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    vtu_write_impl(rPath, rMesh, binary, codec, /*Appended=*/false);
}

void write_vtu_appended(const std::string& rPath, const Mesh& rMesh, detail::VtkCodec codec) {
    vtu_write_impl(rPath, rMesh, /*binary=*/true, codec, /*Appended=*/true);
}

}  // namespace meshioplusplus
