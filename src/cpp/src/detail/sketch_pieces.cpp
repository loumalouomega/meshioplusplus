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
// A fit_curves sketch in world coordinates. See detail/sketch_pieces.hpp.

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "sketch_pieces.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

bool has_sketch(const Mesh& rMesh) {
    return rMesh.HasFieldData("sketch:type") && rMesh.HasFieldData("sketch:params") &&
           rMesh.HasFieldData("sketch:frame") && rMesh.HasFieldData("sketch:frames");
}

std::vector<SketchWorldPiece> sketch_world_pieces(const Mesh& rMesh, const char* pPrefix) {
    std::vector<SketchWorldPiece> out;
    if (!has_sketch(rMesh))
        return out;
    const NDArray& types = rMesh.FieldData("sketch:type");
    const NDArray& params = rMesh.FieldData("sketch:params");
    const NDArray& frame_of = rMesh.FieldData("sketch:frame");
    const NDArray& frames = rMesh.FieldData("sketch:frames");
    if (params.Size() != 10 * types.Size() || frame_of.Size() != types.Size())
        throw WriteError(std::string(pPrefix) + "the sketch:* field data do not match");
    out.reserve(types.Size());
    for (std::size_t k = 0; k < types.Size(); ++k) {
        SketchWorldPiece piece;
        piece.mType = read_int(types, k);
        const std::int64_t f = read_int(frame_of, k);
        if (f < 0 || static_cast<std::size_t>(16 * (f + 1)) > frames.Size())
            throw WriteError(std::string(pPrefix) + "a sketch piece names a missing frame");
        const auto m = [&](std::size_t R, std::size_t C) {
            return read_double(frames, static_cast<std::size_t>(16 * f) + 4 * R + C);
        };
        for (std::size_t r = 0; r < 3; ++r) {
            piece.mE1[r] = m(r, 0);
            piece.mE2[r] = m(r, 1);
            piece.mNormal[r] = m(r, 2);
            piece.mOrigin[r] = m(r, 3);
        }
        double p[10];
        for (std::size_t q = 0; q < 10; ++q)
            p[q] = read_double(params, 10 * k + q);
        const auto world = [&](double U, double V) {
            std::array<double, 3> w{};
            for (std::size_t r = 0; r < 3; ++r)
                w[r] = piece.mOrigin[r] + U * piece.mE1[r] + V * piece.mE2[r];
            return w;
        };
        piece.mStart = world(p[0], p[1]);
        piece.mEnd = world(p[2], p[3]);
        piece.mCentre = world(p[4], p[5]);
        piece.mRadius = p[6];
        piece.mSweepDeg = p[7];
        out.push_back(piece);
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
