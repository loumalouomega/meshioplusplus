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

// Two links: absent affine/empty pairs, then translation and ordered duplicate
// pairs. File node tags (30,10,40,20) map to rows (0,1,2,3).
inline std::string gmsh_periodic_section(int Version, bool Binary, int Width = 8) {
    std::string bytes = "$Periodic\n";
    const bool text = !Binary || Version == 22;
    auto put_i = [&](std::int32_t v) {
        if (text)
            bytes += std::to_string(v) + " ";
        else
            bytes.append(reinterpret_cast<const char*>(&v), 4);
    };
    auto put_size = [&](std::uint64_t v) {
        if (text)
            bytes += std::to_string(v) + "\n";
        else
            bytes.append(reinterpret_cast<const char*>(&v), Width);
    };
    auto put_double = [&](double v) {
        if (text)
            bytes += std::to_string(v) + " ";
        else
            bytes.append(reinterpret_cast<const char*>(&v), 8);
    };
    if (Version == 41)
        put_size(2);
    else {
        put_i(2);
        if (text)
            bytes += "\n";
    }
    for (int i = 0; i < 2; ++i) {
        for (int tag : {1, 11 + i, 33})
            put_i(tag);
        if (text)
            bytes += "\n";
        if (Version == 41)
            put_size(i ? 16 : 0);
        else if (i) {
            if (text)
                bytes += "Affine ";
            else
                put_size(UINT64_MAX);
        }
        if (i) {
            for (int c = 0; c < 16; ++c)
                put_double(c % 5 == 0 || c == 3 ? 1.0 : 0.0);
            if (text)
                bytes += "\n";
        }
        put_size(i ? 3 : 0);
        if (i)
            for (int tag : {30, 10, 40, 20, 30, 10}) {
                if (Version == 41)
                    put_size(tag);
                else
                    put_i(tag);
            }
        if (text)
            bytes += "\n";
    }
    bytes += "\n$EndPeriodic\n";
    return bytes;
}

inline std::string gmsh_periodic_fixture(int Version, bool Binary, int Width = 8) {
    if (Version == 40) {
        auto base = gmsh40_fixture(Binary, Width);
        // The 4.1 writer cannot synthesize one shared physical entity from
        // multiple 4.0 curve blocks. Keep raw physical tags but no named regions
        // here; the separate 4.0 fixture tests their membership.
        const auto start = base.find("$PhysicalNames");
        const auto end = base.find("$EndPhysicalNames\n") + 18;
        base.erase(start, end - start);
        return base + gmsh_periodic_section(40, Binary, Width);
    }
    std::string bytes = "$MeshFormat\n" + std::string(Version == 22 ? "2.2" : "4.1") + " " +
                        (Binary ? "1" : "0") + " " + std::to_string(Width) + "\n";
    auto put_i = [&](std::int32_t v) { bytes.append(reinterpret_cast<const char*>(&v), 4); };
    auto put_size = [&](std::uint64_t v) {
        bytes.append(reinterpret_cast<const char*>(&v), Width);
    };
    auto put_double = [&](double v) { bytes.append(reinterpret_cast<const char*>(&v), 8); };
    if (Binary) {
        put_i(1);
        bytes += "\n";
    }
    bytes += "$EndMeshFormat\n$Nodes\n";
    if (!Binary) {
        if (Version == 22)
            bytes += "4\n30 1 0 0\n10 0 0 0\n40 1 1 0\n20 0 1 0\n";
        else
            bytes += "1 4 10 40\n2 22 0 4\n30 10 40 20\n1 0 0\n0 0 0\n1 1 0\n0 1 0\n";
    } else {
        if (Version == 22)
            bytes += "4\n";
        else {
            for (int v : {1, 4, 10, 40})
                put_size(v);
            for (int v : {2, 22, 0})
                put_i(v);
            put_size(4);
            for (int v : {30, 10, 40, 20})
                put_size(v);
        }
        const int tags[4] = {30, 10, 40, 20};
        const double points[4][3] = {{1, 0, 0}, {0, 0, 0}, {1, 1, 0}, {0, 1, 0}};
        for (int i = 0; i < 4; ++i) {
            if (Version == 22)
                put_i(tags[i]);
            for (double v : points[i])
                put_double(v);
        }
        bytes += "\n";
    }
    bytes += "$EndNodes\n$Elements\n";
    if (!Binary)
        bytes +=
            Version == 22 ? "1\n1 3 2 7 22 10 30 40 20\n" : "1 1 1 1\n2 22 3 1\n1 10 30 40 20\n";
    else {
        if (Version == 22) {
            bytes += "1\n";
            for (int v : {3, 1, 2, 1, 7, 22, 10, 30, 40, 20})
                put_i(v);
        } else {
            for (int v : {1, 1, 1, 1})
                put_size(v);
            for (int v : {2, 22, 3})
                put_i(v);
            put_size(1);
            for (int v : {1, 10, 30, 40, 20})
                put_size(v);
        }
        bytes += "\n";
    }
    bytes += "$EndElements\n";
    return bytes + gmsh_periodic_section(Version, Binary, Width);
}

}  // namespace mt
