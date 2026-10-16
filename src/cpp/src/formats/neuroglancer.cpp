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

// System includes
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/neuroglancer.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "../detail/open_source.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

std::uint32_t ng_le32(const char* pP) {
    const auto* u = reinterpret_cast<const unsigned char*>(pP);
    return static_cast<std::uint32_t>(u[0]) | (static_cast<std::uint32_t>(u[1]) << 8) |
           (static_cast<std::uint32_t>(u[2]) << 16) | (static_cast<std::uint32_t>(u[3]) << 24);
}

void ng_put_le32(std::string& rOut, std::uint32_t V) {
    for (int k = 0; k < 4; ++k)
        rOut.push_back(static_cast<char>((V >> (8 * k)) & 0xffu));
}

}  // namespace

Mesh read_neuroglancer(const std::string& rPath) {
    const detail::FileSource source = detail::open_source(rPath, "Could not open file: " + rPath);
    const std::string_view buf = source.View();
    if (buf.size() < 4)
        throw ReadError("Neuroglancer: the file is shorter than its vertex count");
    const std::size_t n =
        detail::checked_count(ng_le32(buf.data()), (buf.size() - 4) / 12, "Neuroglancer",
                              "vertex");
    const std::size_t tri_bytes = buf.size() - 4 - n * 12;
    if (tri_bytes % 12 != 0)
        throw ReadError("Neuroglancer: the size of the precomputed mesh data is not adequate");

    NDArray points = NDArray::Uninit(DType::Float32, {n, 3});
    float* pp = points.As<float>();
    for (std::size_t i = 0; i < n * 3; ++i) {
        const std::uint32_t bits = ng_le32(buf.data() + 4 + 4 * i);
        std::memcpy(pp + i, &bits, 4);
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(points));

    const std::size_t nt = tri_bytes / 12;
    NDArray cells = NDArray::Uninit(DType::Int64, {nt, 3});
    std::int64_t* cp = cells.As<std::int64_t>();
    const char* tri = buf.data() + 4 + n * 12;
    for (std::size_t i = 0; i < nt * 3; ++i) {
        const std::uint32_t id = ng_le32(tri + 4 * i);
        if (id >= n)
            throw ReadError("Neuroglancer: the mesh references nonexistent vertices");
        cp[i] = id;
    }
    mesh.AddCellBlock("triangle", std::move(cells));
    return mesh;
}

void write_neuroglancer(const std::string& rPath, const Mesh& rMesh) {
    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    if (n > UINT32_MAX)
        throw WriteError("Neuroglancer: too many points for a 32-bit vertex count");

    std::string out;
    ng_put_le32(out, static_cast<std::uint32_t>(n));
    const detail::DoubleView values(points);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            const float f = k < dim ? static_cast<float>(values.Data()[i * dim + k]) : 0.0f;
            std::uint32_t bits;
            std::memcpy(&bits, &f, 4);
            ng_put_le32(out, bits);
        }

    std::size_t skipped = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "triangle") {
            ++skipped;
            continue;
        }
        const NDArray& data = cb.Conn();
        const detail::Int64View ids(data);
        for (std::size_t i = 0; i < data.Size(); ++i) {
            const std::int64_t id = ids.Data()[i];
            if (id < 0 || static_cast<std::uint64_t>(id) > UINT32_MAX)
                throw WriteError("Neuroglancer: a point index does not fit in 32 bits");
            ng_put_le32(out, static_cast<std::uint32_t>(id));
        }
    }
    if (skipped)
        log::warn("Neuroglancer only supports triangle cells; {} other cell block(s) skipped.",
                  skipped);

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
    if (!os)
        throw WriteError("Neuroglancer: could not write " + rPath);
}

}  // namespace meshioplusplus
