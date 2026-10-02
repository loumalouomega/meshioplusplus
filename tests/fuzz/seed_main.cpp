// SPDX-License-Identifier: MIT
// Project-generated positive seeds: no Python, LFS or third-party fixture licence.
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "fuzz_bundle.hpp"
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
                // Companion-file readers need one bundle input, not separate
                // files: pack what the writer left (an XDMF `.xdmf` plus its
                // `.h5`/`.bin`, or a VTX `.bp` directory) so OSS-Fuzz and the
                // campaign start from a valid reconstruction.
                if (format == "xdmf" || format == "vtx") {
                    std::vector<meshioplusplus::FuzzBundleEntry> bundle;
                    for (const auto& entry : std::filesystem::recursive_directory_iterator(dir)) {
                        if (!entry.is_regular_file())
                            continue;
                        const auto rel = std::filesystem::relative(entry.path(), dir);
                        if (std::filesystem::file_size(entry.path()) > 262144)
                            continue;
                        std::ifstream in(entry.path(), std::ios::binary);
                        std::vector<std::uint8_t> data((std::istreambuf_iterator<char>(in)),
                                                       std::istreambuf_iterator<char>());
                        std::string name = rel.string();
                        if (format == "vtx" && name.rfind("seed.bp", 0) == 0)
                            name = "input.bp" + name.substr(7);
                        if (meshioplusplus::fuzz_bundle_name_ok(name))
                            bundle.emplace_back(name, std::move(data));
                        if (bundle.size() >= meshioplusplus::kFuzzBundleMaxEntries)
                            break;
                    }
                    if (!bundle.empty()) {
                        const auto blob = meshioplusplus::fuzz_bundle_encode(bundle);
                        if (blob.size() <= 262144) {
                            const auto out =
                                root / format /
                                ((surface ? "surface-" : "volume-") + std::string("bundle.miob"));
                            std::ofstream f(out, std::ios::binary);
                            f.write(reinterpret_cast<const char*>(blob.data()),
                                    static_cast<std::streamsize>(blob.size()));
                        }
                    }
                }
            } catch (const std::exception& e) {
                std::fprintf(stderr, "seed %s (%s): %s\n", format.c_str(),
                             surface ? "surface" : "volume", e.what());
            }
        }
    }
}
