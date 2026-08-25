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
// Conservative content-based format detection: match only unambiguous leading
// signatures, return "" otherwise. Signatures shared by several formats (the
// generic HDF5 magic, a headerless binary STL) are intentionally NOT claimed.

// System includes
#include <algorithm>
#include <cctype>
#include <fstream>
#include <string>

// Project includes
#include "meshioplusplus/operations/sniff.hpp"

namespace meshioplusplus {

namespace {

// Whether `hay` starts with `needle`.
bool sniff_starts_with(const std::string& rHay, const char* pNeedle) {
    const std::string needle(pNeedle);
    return rHay.size() >= needle.size() && rHay.compare(0, needle.size(), needle) == 0;
}

// Whether `hay` contains `needle`.
bool sniff_contains(const std::string& rHay, const char* pNeedle) {
    return rHay.find(pNeedle) != std::string::npos;
}

// The prefix with leading ASCII whitespace removed.
std::string sniff_lstrip(const std::string& rIn) {
    std::size_t i = 0;
    while (i < rIn.size() && std::isspace(static_cast<unsigned char>(rIn[i])) != 0)
        ++i;
    return rIn.substr(i);
}

}  // namespace

std::string sniff_format(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    if (!in)
        return "";
    char buf[512];
    in.read(buf, sizeof(buf));
    const std::string head(buf, static_cast<std::size_t>(in.gcount()));
    if (head.empty())
        return "";
    const std::string stripped = sniff_lstrip(head);

    // --- binary magics ---
    // VTK XML formats begin (possibly after a BOM/whitespace) with "<?xml" or
    // directly a "<VTKFile" element carrying the grid type.
    if (sniff_contains(head, "VTKFile")) {
        if (sniff_contains(head, "UnstructuredGrid"))
            return "vtu";
        if (sniff_contains(head, "PolyData"))
            return "vtp";
        // Checked last of the three: the grid-type strings are disjoint, but a
        // future dataset type could contain another as a substring, and the
        // cheapest defence is to keep the most recently added one from
        // shadowing anything.
        if (sniff_contains(head, "ImageData"))
            return "vti";
    }
    if (sniff_starts_with(stripped, "<Xdmf") || sniff_contains(head, "<Xdmf"))
        return "xdmf";

    // --- line/text signatures ---
    if (sniff_starts_with(stripped, "# vtk DataFile"))
        return "vtk";
    if (sniff_starts_with(stripped, "$MeshFormat"))
        return "gmsh";
    // GiD postprocess. The results file is self-identifying. The geometry file
    // is matched on `MESH "` -- keyword, space AND opening quote -- because a
    // bare "MESH " prefix is exactly the generic English token this file's own
    // contract warns against claiming; the quote is what gidpost always writes
    // and what makes the match unambiguous. Neither `.post.bin` (a deflated
    // stream, no stable leading signature) nor `.post.h5` (the generic HDF5
    // magic, which this file deliberately never claims) is sniffable.
    if (sniff_starts_with(stripped, "GiD Post Results File"))
        return "gid";
    if (sniff_starts_with(stripped, "MESH \""))
        return "gid";
    // PLY: "ply" on its own first line.
    if (sniff_starts_with(stripped, "ply\n") || sniff_starts_with(stripped, "ply\r") ||
        stripped == "ply")
        return "ply";
    // OFF variants (OFF / COFF / NOFF / STOFF ...): a token ending in "OFF".
    if (sniff_starts_with(stripped, "OFF") || sniff_starts_with(stripped, "COFF") ||
        sniff_starts_with(stripped, "NOFF") || sniff_starts_with(stripped, "STOFF"))
        return "off";
    // ASCII STL.
    if (sniff_starts_with(stripped, "solid "))
        return "stl";
    // Abaqus input decks start with a keyword line "*Heading"/"*Node"/"*NODE".
    {
        std::string upper = stripped.substr(0, 8);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (sniff_starts_with(upper, "*HEADING") || sniff_starts_with(upper, "*NODE"))
            return "abaqus";
    }

    return "";
}

}  // namespace meshioplusplus
