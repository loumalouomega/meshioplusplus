// SPDX-License-Identifier: MIT
#pragma once

#include <cstdint>
#include <string>

namespace mt {

// Hand-built 4.0 records, independent of any writer. Sparse/out-of-order
// node tags, distinct same-type blocks, physical names and a point entity
// (4.0 stores six bounding-box doubles even for points).
inline std::string gmsh40_fixture(bool Binary, int CountBytes = 8) {
    if (!Binary)
        return R"($MeshFormat
4.0 0 8
$EndMeshFormat
$PhysicalNames
2
1 8 "edge"
2 7 "plate"
$EndPhysicalNames
$Entities
1 2 1 0
5 0 0 0 0 0 0 0
11 0 0 0 1 0 0 1 8 2 5 -6
33 0 0 0 0 1 0 1 8 0
22 0 0 0 1 1 0 1 7 2 11 -33
$EndEntities
$Nodes
2 4
11 1 0 2
30 1 0 0
10 0 0 0
22 2 0 2
40 1 1 0
20 0 1 0
$EndNodes
$Elements
3 4
11 1 1 1
90 10 30
22 2 2 2
80 10 30 40
60 10 40 20
33 1 1 1
50 20 10
$EndElements
)";

    std::string bytes = "$MeshFormat\n4.0 1 8\n";
    auto put_i = [&](std::int32_t V) { bytes.append(reinterpret_cast<const char*>(&V), 4); };
    auto put_count = [&](std::uint64_t V) {
        if (CountBytes == 4) {
            const auto v = static_cast<std::uint32_t>(V);
            bytes.append(reinterpret_cast<const char*>(&v), 4);
        } else {
            bytes.append(reinterpret_cast<const char*>(&V), 8);
        }
    };
    auto put_d = [&](double V) { bytes.append(reinterpret_cast<const char*>(&V), 8); };
    put_i(1);
    bytes +=
        "\n$EndMeshFormat\n$PhysicalNames\n2\n1 8 \"edge\"\n2 7 \"plate\"\n$EndPhysicalNames\n";
    bytes += "$Entities\n";
    for (int v : {1, 2, 1, 0})
        put_count(v);
    for (int tag : {5, 11, 33, 22}) {
        put_i(tag);
        for (int c = 0; c < 6; ++c)
            put_d(0.0);
        put_count(tag == 5 ? 0 : 1);
        if (tag != 5)
            put_i(tag == 22 ? 7 : 8);
        if (tag != 5) {
            put_count(tag == 33 ? 0 : 2);
            if (tag != 33) {
                put_i(tag == 11 ? 5 : 11);
                put_i(tag == 11 ? -6 : -33);
            }
        }
    }
    bytes += "\n$EndEntities\n$Nodes\n";
    put_count(2);
    put_count(4);
    const int tags[4] = {30, 10, 40, 20};
    const double points[4][3] = {{1, 0, 0}, {0, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    for (int block = 0; block < 2; ++block) {
        put_i(block == 0 ? 11 : 22);
        put_i(block + 1);
        put_i(0);
        put_count(2);
        for (int i = block * 2; i < block * 2 + 2; ++i) {
            put_i(tags[i]);
            for (double v : points[i])
                put_d(v);
        }
    }
    bytes += "\n$EndNodes\n$Elements\n";
    put_count(3);
    put_count(4);
    const int headers[3][4] = {{11, 1, 1, 1}, {22, 2, 2, 2}, {33, 1, 1, 1}};
    const int rows[4][4] = {{90, 10, 30, 0}, {80, 10, 30, 40}, {60, 10, 40, 20}, {50, 20, 10, 0}};
    int row = 0;
    for (const auto& header : headers) {
        for (int c = 0; c < 3; ++c)
            put_i(header[c]);
        put_count(header[3]);
        for (int e = 0; e < header[3]; ++e, ++row)
            for (int c = 0; c < (header[2] == 1 ? 3 : 4); ++c)
                put_i(rows[row][c]);
    }
    bytes += "\n$EndElements\n";
    return bytes;
}

}  // namespace mt
