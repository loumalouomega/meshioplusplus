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
 * @file mesh_digest.hpp
 * @brief A byte-level digest of a mesh, for determinism checks.
 *
 * Private to the benchmark and the C++ tests: not installed, not part of the
 * API or the ABI.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"

namespace meshioplusplus::bench {

/**
 * @brief 64-bit FNV-1a over everything a mesh (or an operation result)
 * carries, in a fixed order: points, every block's type and connectivity
 * (ragged rows and polyhedron faces included), point/cell/field data in
 * sorted-name order, and regions. Feed index maps and counters through
 * `Array`/`Arrays`/`U64` after `Of`.
 *
 * Two results with equal digests are byte-identical for every practical
 * purpose; this is what `meshioplusplus_bench_ops --hash` prints and what the
 * C++ determinism tests compare (roadmap §3).
 */
class MeshDigest {
public:
    void Bytes(const void* pData, std::size_t n) {
        const auto* b = static_cast<const unsigned char*>(pData);
        for (std::size_t i = 0; i < n; ++i) {
            mH ^= b[i];
            mH *= 1099511628211ull;
        }
    }
    void U64(std::uint64_t v) { Bytes(&v, sizeof v); }
    void Str(const std::string& rS) {
        U64(rS.size());
        Bytes(rS.data(), rS.size());
    }
    void Array(const meshioplusplus::NDArray& rA) {
        U64(static_cast<std::uint64_t>(rA.Dtype()));
        U64(rA.Shape().size());
        for (std::size_t d : rA.Shape())
            U64(d);
        Bytes(rA.Data(), rA.Nbytes());
    }
    void Arrays(const std::vector<meshioplusplus::NDArray>& rAs) {
        U64(rAs.size());
        for (const meshioplusplus::NDArray& a : rAs)
            Array(a);
    }
    void Of(const meshioplusplus::Mesh& rM) {
        Array(rM.Points());
        U64(rM.NumCellBlocks());
        for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b) {
            const auto cb = rM.Cells(b);
            Str(cb.Type());
            U64(cb.NumCells());
            if (cb.IsPolyhedron()) {
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    U64(cb.NumFaces(c));
                    for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                        const auto face = cb.Face(c, f);
                        U64(face.second);
                        Bytes(face.first, face.second * sizeof(std::int64_t));
                    }
                }
            } else if (cb.IsRagged()) {
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    U64(cb.RowSize(c));
                    Bytes(cb.Row(c), cb.RowSize(c) * sizeof(std::int64_t));
                }
            } else {
                Array(cb.Conn());
            }
        }
        for (const std::string& n : rM.PointDataNames()) {
            Str(n);
            Array(rM.PointData(n));
        }
        for (const std::string& n : rM.CellDataNames()) {
            Str(n);
            for (std::size_t b = 0; b < rM.NumCellBlocks(); ++b)
                Array(rM.CellData(n, b));
        }
        for (const std::string& n : rM.FieldDataNames()) {
            Str(n);
            Array(rM.FieldData(n));
        }
        U64(rM.NumRegions());
        for (std::size_t r = 0; r < rM.NumRegions(); ++r) {
            const meshioplusplus::Region& rg = rM.Region(r);
            Str(rg.mName);
            U64(static_cast<std::uint64_t>(rg.mKind));
            U64(static_cast<std::uint64_t>(static_cast<std::int64_t>(rg.mDim)));
            U64(static_cast<std::uint64_t>(static_cast<std::int64_t>(rg.mTag)));
            Array(rg.mEntries);
        }
    }
    std::uint64_t Value() const { return mH; }

private:
    std::uint64_t mH = 14695981039346656037ull;
};

}  // namespace meshioplusplus::bench
