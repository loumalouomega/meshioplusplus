// SPDX-License-Identifier: MIT
// Project-generated positive seeds: no Python, LFS or third-party fixture licence.
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "fuzz_bundle.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

// The seed kinds: a tetra, a triangle, and -- only for a format whose writer took
// neither (a drawing such as DXF) -- a polyline; and a writer that takes none of
// them (PLOT3D's lattice blocks) -- one 2x2x2 lattice block.
enum class SeedKind { Volume, Surface, Line, Lattice };

const char* seed_label(SeedKind Kind) {
    return Kind == SeedKind::Volume    ? "volume"
           : Kind == SeedKind::Surface ? "surface"
           : Kind == SeedKind::Line    ? "line"
                                       : "lattice";
}

meshioplusplus::Mesh seed_mesh(SeedKind Kind) {
    using namespace meshioplusplus;
    Mesh mesh;
    if (Kind == SeedKind::Lattice) {
        // One 2x2x2 block for the writers that need a lattice (PLOT3D's plot3d:dims).
        NDArray points = NDArray::Uninit(DType::Float64, {8, 3});
        for (std::size_t p = 0; p < 8; ++p)
            for (std::size_t c = 0; c < 3; ++c)
                points.As<double>()[p * 3 + c] = static_cast<double>((p >> c) & 1u);
        mesh.AssignPoints(std::move(points));
        NDArray cells = NDArray::Uninit(DType::Int64, {1, 8});
        const std::int64_t conn[] = {0, 1, 3, 2, 4, 5, 7, 6};
        std::copy(std::begin(conn), std::end(conn), cells.As<std::int64_t>());
        mesh.AddCellBlock("hexahedron", std::move(cells));
        NDArray dims = NDArray::Uninit(DType::Int64, {1, 3});
        std::fill_n(dims.As<std::int64_t>(), 3, std::int64_t{2});
        mesh.AddFieldData("plot3d:dims", std::move(dims));
        return mesh;
    }
    NDArray points = NDArray::Uninit(DType::Float64, {4, 3});
    const double coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::copy(std::begin(coords), std::end(coords), points.As<double>());
    mesh.AssignPoints(std::move(points));
    if (Kind == SeedKind::Line) {
        NDArray cells = NDArray::Uninit(DType::Int64, {3, 2});
        const std::int64_t conn[] = {0, 1, 1, 2, 2, 3};
        std::copy(std::begin(conn), std::end(conn), cells.As<std::int64_t>());
        mesh.AddCellBlock("line", std::move(cells));
        return mesh;
    }
    const bool surface = Kind == SeedKind::Surface;
    NDArray cells = NDArray::Uninit(DType::Int64, {1, surface ? 3u : 4u});
    for (std::size_t i = 0; i < cells.Size(); ++i)
        cells.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellBlock(surface ? "triangle" : "tetra", std::move(cells));
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
        bool seeded = false;
        for (SeedKind kind :
             {SeedKind::Volume, SeedKind::Surface, SeedKind::Line, SeedKind::Lattice}) {
            if ((kind == SeedKind::Line || kind == SeedKind::Lattice) && seeded)
                continue;
            const std::string label = seed_label(kind);
            const auto dir = root / format / label;
            std::filesystem::create_directories(dir);
            const auto path = dir / ("seed" + extension);
            try {
                writer(path.string(), seed_mesh(kind));
                // Validate through the native reader before supplying a seed.
                (void)meshioplusplus::registry_readers().at(format)(path.string());
                seeded = true;
                for (const auto& entry : std::filesystem::directory_iterator(dir)) {
                    if (entry.is_regular_file() &&
                        std::filesystem::file_size(entry.path()) <= 262144)
                        std::filesystem::copy_file(
                            entry.path(),
                            root / format / (label + "-" + entry.path().filename().string()),
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
                            const auto out = root / format / (label + "-bundle.miob");
                            std::ofstream f(out, std::ios::binary);
                            f.write(reinterpret_cast<const char*>(blob.data()),
                                    static_cast<std::streamsize>(blob.size()));
                        }
                    }
                }
            } catch (const std::exception& e) {
                std::fprintf(stderr, "seed %s (%s): %s\n", format.c_str(), label.c_str(), e.what());
            }
        }
    }
}
