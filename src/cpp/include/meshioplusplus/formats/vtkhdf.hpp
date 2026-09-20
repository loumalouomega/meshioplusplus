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
 * @file vtkhdf.hpp
 * @brief VTKHDF (`.vtkhdf`) C++ reader/writer: Kitware's HDF5-based VTK format,
 *        the one ParaView-native format that carries time, partitions and fields
 *        in a single file.
 *
 * Reads and writes `Type` `UnstructuredGrid`, `PolyData`,
 * `PartitionedDataSetCollection` and `MultiBlockDataSet`. `ImageData`,
 * `OverlappingAMR`, `HyperTreeGrid`, `Table`, `RectilinearGrid` and
 * `StructuredGrid` are refused by name. Everything below the API is specified by
 * `doc/formats/vtkhdf.md`; the layout facts it depends on were measured against
 * VTK 9.7's own `vtkHDFWriter`/`vtkHDFReader`, not taken from the prose spec.
 *
 * ### Reading
 *
 * - A file with several partitions is **merged** into one mesh with one
 *   `RegionKind::Cell` region per piece (`piece_<i>`; a composite's blocks are
 *   named by the block). `ReadOptions::mPieceSet` / `mPiece` selects one piece
 *   instead, and then no regions are attached.
 * - `ReadOptions::mTimeStep` selects a step of a transient file; the step's time
 *   is attached as `field_data["meshio:time"]`.
 * - Polyhedra are bucketed into `polyhedron<N>` blocks (N = unique node count)
 *   **per piece**, so a piece's cells are always one contiguous range and reading
 *   piece `k` alone reproduces region `k` of the merged read.
 * - The reader dispatches on which groups exist, never on the declared minor
 *   version: files in the wild declare sloppy versions.
 * - A transient *composite* (a `Steps` group with more than one step on a
 *   `PartitionedDataSetCollection`/`MultiBlockDataSet`) is refused by name.
 *
 * ### Writing
 *
 * The version written is the oldest that covers the features used (2.0 for a plain
 * mesh, 2.1 for a composite, 2.5 with polyhedra) unless pinned. `Type` is a
 * fixed-length ASCII attribute. `PolyData` is never chosen automatically and
 * refuses any cell that is not a vertex, line, triangle, quad or polygon.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief The VTKHDF dataset kind written by `write_vtkhdf`. */
enum class VtkhdfType {
    UnstructuredGrid,
    PolyData,
    PartitionedDataSetCollection,
    MultiBlockDataSet,
};

/** @brief A VTKHDF version to declare; `{0, 0}` (the default) means "automatic". */
struct VtkhdfVersion {
    int mMajor = 0;
    int mMinor = 0;
};

/**
 * @brief Write a static VTKHDF file.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param GzipLevel deflate level 0-9 for datasets of at least 4 KiB, or negative
 *        for none; smaller datasets are never compressed, because chunk metadata
 *        would outweigh any saving
 * @param Type the dataset kind. The two composite kinds write one block per cell
 *        region (when the regions partition the cells) or per cell block, and
 *        replicate `field_data` onto every block (`vtkHDFReader` fails on any other
 *        root group)
 * @param Version the version to declare; automatic when `{0, 0}`. A version too
 *        old for the content is a `WriteError` naming the feature.
 * @throws WriteError for an unwritable mesh (a cell type with no VTK id, a
 *         non-numeric array, a PolyData volume cell, an empty composite) or path.
 */
MESHIOPLUSPLUS_API void write_vtkhdf(const std::string& rPath, const Mesh& rMesh, int GzipLevel = 4,
                                     VtkhdfType Type = VtkhdfType::UnstructuredGrid,
                                     VtkhdfVersion Version = {});

/** @brief Read a VTKHDF file with default options: every piece merged, step 0. */
MESHIOPLUSPLUS_API Mesh read_vtkhdf(const std::string& rPath);

/**
 * @brief Read a VTKHDF file honouring `rOpts`.
 *
 * Honoured: `mPointsOnly`, `mDataArrays`, `mTimeStep`, `mPiece`/`mPieceSet` and
 * `mLenient` (skips cells with no meshio++ type -- poly-vertex, poly-line,
 * triangle strips -- with a warning instead of throwing).
 * @throws ReadError if the file is not VTKHDF, is an unsupported type or version,
 *         or a step/piece is out of range.
 */
MESHIOPLUSPLUS_API Mesh read_vtkhdf(const std::string& rPath, const ReadOptions& rOpts);

/**
 * @brief Summarize a VTKHDF file without loading its heavy arrays.
 *
 * Reads only `Steps/Values`, the `NumberOf*` counters, `Types` and the data
 * names, so `mFellBackToFullRead` is false for a polyhedron-free
 * `UnstructuredGrid`. `PolyData`, composites and polyhedral files fall back to a
 * full read (their block structure needs the arrays). `mTimeValues` is filled in
 * either case.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vtkhdf_metadata(const std::string& rPath,
                                                     const ReadOptions& rOpts);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
