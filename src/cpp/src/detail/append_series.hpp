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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file detail/append_series.hpp
 * @brief Shared step-append helper for single-mesh series writers that write
 * step 0 in full and append only data sections extracted from a scratch write
 * for subsequent steps (Gmsh, UNV; roadmap 6.4).
 *
 * Core-private (not installed).
 */

// System includes
#include <cstddef>
#include <cstdio>
#include <ios>
#include <iterator>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::detail {

template <class WriteFn, class ExtractFn>
void append_series_write_step(const std::string& rPath, const Mesh& rMesh, std::size_t& rSteps,
                              std::size_t& rPoints, std::vector<std::size_t>& rBlockCells,
                              const char* pSeriesLabel, WriteFn&& rWriteFull,
                              ExtractFn&& rExtractData) {
    if (rSteps == 0) {
        rWriteFull(rPath);
        rPoints = rMesh.NumPoints();
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            rBlockCells.push_back(rMesh.Cells(b).NumCells());
    } else {
        bool same = rMesh.NumPoints() == rPoints && rMesh.NumCellBlocks() == rBlockCells.size();
        for (std::size_t b = 0; same && b < rBlockCells.size(); ++b)
            same = rMesh.Cells(b).NumCells() == rBlockCells[b];
        if (!same)
            throw WriteError(std::string(pSeriesLabel) + ": step " + std::to_string(rSteps) +
                             "'s points or cells differ from the first step's; a file holds one "
                             "mesh -- write one file per step with '{step}'");
        const std::string scratch = rPath + ".step.tmp";
        rWriteFull(scratch);
        std::string bytes;
        {
            auto in = make_classic_ifstream(scratch, std::ios::binary);
            bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
        }
        std::remove(scratch.c_str());
        const std::string extra = rExtractData(bytes);
        auto out = make_classic_ofstream(rPath, std::ios::binary | std::ios::app);
        if (!out)
            throw WriteError("Could not open file for writing: " + rPath);
        out.write(extra.data(), static_cast<std::streamsize>(extra.size()));
    }
    ++rSteps;
}

}  // namespace meshioplusplus::detail
