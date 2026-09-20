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
 * @file vtk_xml.hpp
 * @brief Shared VTK-XML `<DataArray>` helpers used by the VTU
 * (UnstructuredGrid) and VTP (PolyData) readers/writers.
 *
 * Both formats are the same XML container — `<VTKFile>` with per-array
 * `<DataArray>` elements in ascii or base64 "binary" encoding (raw or
 * zlib-compressed, framed by `detail/vtu_binary.hpp`) — differing only in
 * the grid element in between. This header holds the container-level pieces:
 * the DType <-> VTK type-name mapping (`vtu_type_str`/`dtype_from_vtu`),
 * ASCII float/array emission (`vtu_ascii_double`/`vtu_ascii_ndarray`),
 * DataArray text parsing for both encodings (`vtu_parse_ascii` /
 * `vtu_parse_binary`), and the `vtu_to_int64` widening helper.
 *
 * Deliberately pugixml-free (the single-header amalgamation only bundles
 * pugixml in its implementation section): the readers each keep a thin local
 * wrapper that pulls the `format`/`type`/`NumberOfComponents` attributes off
 * a `<DataArray>` node and dispatches to `vtu_parse_ascii`/`vtu_parse_binary`.
 *
 * Every function here is called once per data array or once per parsed
 * value inside the ASCII path (where `snprintf`/`strtod` already dominate the
 * cost far past any function-call overhead), so bodies live in
 * `src/cpp/src/detail/vtk_xml.cpp` rather than inline here.
 */

// System includes
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Map an `NDArray` dtype to its VTK-XML type-name string.
 * @param dt The dtype.
 * @return The VTK type name (e.g. `"Float64"`, `"Int32"`).
 */
MESHIOPLUSPLUS_API const char* vtu_type_str(DType dt);

/**
 * @brief Map a VTK-XML type-name string to the `NDArray` dtype.
 * @param rS The VTK type name.
 * @return The matching dtype.
 * @throws ReadError on an unknown type name.
 */
MESHIOPLUSPLUS_API DType dtype_from_vtu(const std::string& rS);

/**
 * @brief Emit one float in VTK's `%.11e` ASCII format followed by a newline.
 * @param rOs Output stream.
 * @param v The value.
 */
MESHIOPLUSPLUS_API void vtu_ascii_double(std::ostream& rOs, double v);

/**
 * @brief Emit a whole `NDArray` in ASCII, one value per line (floats via
 * `vtu_ascii_double`, integers as plain decimals).
 * @param rOs Output stream.
 * @param rA The array.
 */
MESHIOPLUSPLUS_API void vtu_ascii_ndarray(std::ostream& rOs, const NDArray& rA);

/**
 * @brief Store one parsed value into a dtype-erased array slot.
 * @param rA Destination array.
 * @param i Flat index.
 * @param d The value when `rA` is a float dtype.
 * @param v The value when `rA` is an integer dtype.
 */
MESHIOPLUSPLUS_API void vtu_store(NDArray& rA, std::size_t i, double d, std::int64_t v);

/**
 * @brief Parse whitespace-separated ASCII DataArray text into a flat array.
 * @param pText The element text (may be null).
 * @param dt Target dtype (drives float vs integer parsing).
 * @return A 1-D owning array of every parsed value.
 */
MESHIOPLUSPLUS_API NDArray vtu_parse_ascii(const char* pText, DType dt);

/**
 * @brief Trim leading/trailing whitespace from a C string.
 * @param pS The string (may be null).
 * @return The trimmed copy.
 */
MESHIOPLUSPLUS_API std::string vtu_strip(const char* pS);

/**
 * @brief Decode a base64 "binary" DataArray payload into a flat array.
 * @param rText The stripped base64 text.
 * @param dt Target dtype.
 * @param codec block-compression codec recorded in the file.
 * @param hsz Header integer size in bytes (4 for UInt32, 8 for UInt64).
 * @return A 1-D owning array over the decoded bytes.
 */
MESHIOPLUSPLUS_API NDArray vtu_parse_binary(const std::string& rText, DType dt, VtkCodec codec, std::size_t hsz);

/**
 * @brief Widen a dtype-erased integer array to a `std::int64_t` vector.
 * @param rA The array.
 * @return The widened values.
 */
MESHIOPLUSPLUS_API std::vector<std::int64_t> vtu_to_int64(const NDArray& rA);

/**
 * @brief The array to write under @p rName: @p rArray, except that `vtkGhostType` is
 * always `UInt8` on disk.
 *
 * `vtkGhostType` is VTK's reserved ghost-flag name, and a reader only recognises it
 * as the ghost array when it is an unsigned char array (ParaView ignores an `Int64`
 * one entirely). The NATIVE and KRATOS mesh backends hold every integer array as
 * `Int64`, so a ghost array they read or derive would otherwise be written -- and
 * declared by a `.pvtu` index -- as `Int64`. Any other name, and an array that is
 * already `UInt8`, is returned as is; otherwise the values are copied into
 * @p rScratch (each cast to a byte, as a ghost flag is a bit set) and that is returned.
 *
 * @param rName the array's name.
 * @param rArray the array as the mesh holds it.
 * @param rScratch storage for the converted copy; untouched when none is needed.
 */
MESHIOPLUSPLUS_API const NDArray& vtu_disk_array(const std::string& rName, const NDArray& rArray,
                                                 NDArray& rScratch);

/// The dtype `vtu_disk_array` writes @p rName as: `UInt8` for `vtkGhostType`, else @p Dt.
MESHIOPLUSPLUS_API DType vtu_disk_dtype(const std::string& rName, DType Dt);

/**
 * @brief Write one `<FieldData>` array as a complete `<DataArray>` element.
 *
 * A field-data array is one value (or row) per *tuple* of the dataset, so unlike
 * a point or cell array it carries an explicit `NumberOfTuples` (VTK's readers
 * require it there): the leading extent, or 1 for a rank-0 scalar. A rank of two
 * or more also gets `NumberOfComponents`, the product of the trailing extents, so
 * `(n, m)` reads back as `(n, m)`; a higher rank is flattened to that shape.
 * The caller writes the surrounding `<FieldData>` element.
 *
 * @param rOs the output stream.
 * @param rName the array's `Name` attribute (written verbatim, like point data's).
 * @param rArray the array.
 * @param Binary base64-encode the body instead of writing text.
 * @param Codec the block compressor for a binary body; ignored when @p Binary is false.
 */
MESHIOPLUSPLUS_API void vtu_write_field_array(std::ostream& rOs, const std::string& rName,
                                              const NDArray& rArray, bool Binary, VtkCodec Codec);

}  // namespace detail
}  // namespace meshioplusplus
