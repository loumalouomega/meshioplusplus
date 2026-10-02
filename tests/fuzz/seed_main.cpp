// SPDX-License-Identifier: MIT
// Project-generated positive seeds: no Python, LFS or third-party fixture licence.
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

meshioplusplus::Mesh seed_mesh(bool Surface) {
    using namespace meshioplusplus;
    Mesh mesh;
    NDArray points = NDArray::Uninit(DType::Float64, {4, 3});
    const double coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::copy(std::begin(coords), std::end(coords), points.As<double>());
    mesh.AssignPoints(std::move(points));
    NDArray cells = NDArray::Uninit(DType::Int64, {1, Surface ? 3u : 4u});
    for (std::size_t i = 0; i < cells.Size(); ++i)
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellBlock(Surface ? "triangle" : "tetra", std::move(cells));
    NDArray u = NDArray::Uninit(DType::Float64, {4});
    std::fill_n(u.As<double>(), 4, 1.0);
    mesh.AddPointData("u", std::move(u));
    return mesh;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s OUT_DIR\n", argv[0]);
        return 2;
    }
    const std::filesystem::path root(argv[1]);
    for (const auto& [format, writer] : meshioplusplus::registry_writers()) {
        if (!meshioplusplus::registry_readers().count(format))
            continue;
        std::string extension = ".dat";
        for (const auto& [ext, fmt] : meshioplusplus::registry_extension_defaults()) {
            if (fmt == format) {
                extension = ext;
                break;
            }
        }
        for (bool surface : {false, true}) {
            const auto dir = root / format / (surface ? "surface" : "volume");
            std::filesystem::create_directories(dir);
            const auto path = dir / ("seed" + extension);
            try {
                writer(path.string(), seed_mesh(surface));
                // Validate through the native reader before supplying a seed.
                (void)meshioplusplus::registry_readers().at(format)(path.string());
                for (const auto& entry : std::filesystem::directory_iterator(dir)) {
                    if (entry.is_regular_file() &&
                        std::filesystem::file_size(entry.path()) <= 262144)
                        std::filesystem::copy_file(
                            entry.path(),
                            root / format /
                                ((surface ? "surface-" : "volume-") +
                                 entry.path().filename().string()),
                            std::filesystem::copy_options::overwrite_existing);
                }
            } catch (const std::exception& e) {
                std::fprintf(stderr, "seed %s (%s): %s\n", format.c_str(),
                             surface ? "surface" : "volume", e.what());
            }
        }
    }
}
