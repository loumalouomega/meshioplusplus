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
 * @file xdmf_common.hpp
 * @brief XDMF cell-type-name maps and cell-data raw<->blocks conversion,
 * shared between the XDMF format implementation and HMF (which reuses
 * XDMF's topology names and raw cell-data layout).
 *
 * Ported from `src/meshio/xdmf/common.py` and the `raw_from_cell_data` /
 * `cell_data_from_raw` helpers in `src/meshio/_common.py`. XDMF (and HMF)
 * store per-cell-type-block data as one array *per cell type name string* in
 * the XML/XDMF Topology, and store cell_data for a mixed mesh as one
 * concatenated raw array per data name (all cell blocks laid end-to-end)
 * rather than one array per block — `concat_cell_data`/`split_raw_cell_data`
 * are what let this header's callers go between meshio's per-block
 * `cell_data` representation and that concatenated-raw representation.
 *
 * It also owns the pieces of the XDMF *write* path that the single-step writer
 * (`write_xdmf`) and the transient writer (`XdmfTimeSeriesWriter`) both need:
 * the dtype/attribute-type/mixed-index lookups and `DataItemStore`, which
 * stores one array in the chosen heavy-data container and hands back the text a
 * `<DataItem>` element carries. Keeping the payload side here rather than
 * transcribing it means the two writers cannot drift on dataset naming (`data0`,
 * `data1`, ... in one flat namespace), on the sibling `.h5`/`.bin` file naming,
 * or on the inline-XML number formatting. The XML *element* side stays with each
 * writer, because pugixml is a build-only vendored dependency that no installed
 * header may include -- which is also why `DataItemStore` hides its HDF5 handle
 * behind a pimpl rather than naming `h5::Hid` here.
 *
 * Each function here is called once per cell-type-name lookup or once per
 * data array (not per element), so bodies live in
 * `src/cpp/src/detail/xdmf_common.cpp` rather than inline here.
 */

// System includes
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace xdmfcommon {

/**
 * @brief Maps a meshio cell-type name to its XDMF topology type name.
 * @param t meshio cell-type name (e.g. `"triangle"`, `"tetra10"`).
 * @return The corresponding XDMF `TopologyType` string (e.g. `"Triangle"`).
 * @throws WriteError if `t` has no XDMF equivalent.
 */
MESHIOPLUSPLUS_API const char* meshio_to_xdmf(const std::string& rT);

/**
 * @brief Maps an XDMF topology type name to a meshio cell-type name.
 *
 * Accepts both the canonical XDMF spelling and common abbreviations some
 * writers emit (e.g. both `"Hexahedron_20"` and `"Hex_20"` map to
 * `"hexahedron20"`).
 * @param t XDMF `TopologyType` string as found in the file.
 * @return The corresponding meshio cell-type name.
 * @throws ReadError if `t` is not a recognized topology type.
 */
MESHIOPLUSPLUS_API std::string xdmf_to_meshio(const std::string& rT);

/**
 * @brief Concatenates one cell-data name's per-block arrays along axis 0
 * into a single raw array, matching Python's `raw_from_cell_data`.
 *
 * Used when writing XDMF/HMF cell data for a mixed-cell-type mesh: XDMF
 * stores cell data as one flat array per data name (all blocks' rows
 * back-to-back) rather than one array per block.
 * @param rMesh The mesh whose cell data to concatenate.
 * @param rName The cell-data name; must have at least one block, all blocks
 *              sharing dtype/trailing shape.
 * @return A new array with the same trailing shape as the first block and
 *         first dimension equal to the sum of each block's row count.
 */
MESHIOPLUSPLUS_API NDArray concat_cell_data(const Mesh& rMesh, const std::string& rName);

/**
 * @brief Splits a raw, whole-mesh cell-data array (as read from XDMF/HMF)
 * back into one `NDArray` per cell block, matching Python's
 * `cell_data_from_raw`.
 *
 * Inverse of `concat_cell_data`.
 * @param raw The concatenated array covering every cell block's rows,
 *            in cell-block order.
 * @param sizes Row count of each cell block, in the same order the blocks
 *              appear in `raw`; must sum to `raw`'s row count.
 * @return One `NDArray` per entry in `sizes`, each holding that block's slice.
 */
MESHIOPLUSPLUS_API std::vector<NDArray> split_raw_cell_data(const NDArray& rRaw,
                                                            const std::vector<std::size_t>& rSizes);

// ---- write-path helpers (shared by write_xdmf and XdmfTimeSeriesWriter) ----

/**
 * @brief Maps a dtype to the `DataType`/`Precision` attribute pair XDMF spells
 * it with.
 * @param dt The dtype to describe.
 * @return `{"Int"|"UInt"|"Float", "1"|"2"|"4"|"8"}`.
 */
MESHIOPLUSPLUS_API std::pair<const char*, const char*> numpy_to_xdmf_dtype(DType dt);

/**
 * @brief Classifies a data array's shape as an XDMF `AttributeType`.
 * @param rShape The array's shape.
 * @return `"Scalar"`, `"Vector"`, `"Tensor"`, `"Tensor6"` or `"Matrix"`.
 */
MESHIOPLUSPLUS_API std::string attribute_type(const std::vector<std::size_t>& rShape);

/**
 * @brief Maps a meshio cell-type name to its numeric index in a `Mixed`
 * topology array.
 * @param rT meshio cell-type name.
 * @return The XDMF type index (e.g. `0x6` for `"tetra"`).
 * @throws WriteError if the type has no `Mixed` encoding.
 */
MESHIOPLUSPLUS_API int meshio_to_xdmf_index(const std::string& rT);

/**
 * @brief Renders an array's shape as a `Dimensions` attribute value.
 * @param rArr The array.
 * @return Space-separated extents, e.g. `"12 3"`.
 */
MESHIOPLUSPLUS_API std::string dims_string(const NDArray& rArr);

/**
 * @brief Packs every cell block into one flat `Mixed` topology array.
 *
 * Each cell contributes its XDMF type index followed by its node ids;
 * `vertex`/`Polyvertex` and `line`/`Polyline` additionally repeat the index as
 * the point count XDMF requires there (which, for a `line`, is the constant 2
 * the reader insists on).
 * @param rMesh The mesh whose blocks to pack.
 * @param rTotalCells Out: the total cell count, for `NumberOfElements`.
 * @return An `Int64` 1-D array holding the packed topology.
 * @throws WriteError if a block's cell type has no `Mixed` encoding.
 */
MESHIOPLUSPLUS_API NDArray pack_mixed_topology(const Mesh& rMesh, std::size_t& rTotalCells);

/**
 * @brief Stores XDMF heavy data and reports how a `<DataItem>` should reference
 * it.
 *
 * One instance owns one output "container": for `"HDF"` a single sibling
 * `<base>.h5` created lazily on the first store and closed on destruction, for
 * `"Binary"` a series of `<base><n>.bin` files, and for `"XML"` nothing at all
 * (the payload is the returned text). Arrays are numbered `data0`, `data1`, ...
 * across the whole instance's lifetime, so a transient writer's steps all land
 * in one flat namespace exactly as the Python `TimeSeriesWriter` does.
 *
 * `"HDF"` on a build without HDF5 support throws from `Store` naming the CMake
 * option rather than the symbol being absent -- the compiled-out contract every
 * optional-dependency path here follows.
 */
class MESHIOPLUSPLUS_API DataItemStore {
public:
    /**
     * @brief Construct a store.
     * @param rDataFormat `"XML"`, `"Binary"` or `"HDF"`.
     * @param rBase Output path *without* its extension; sibling `.h5`/`.bin`
     *        files are derived from it.
     * @param GzipLevel gzip level for `"HDF"` datasets, negative for none.
     */
    DataItemStore(const std::string& rDataFormat, const std::string& rBase, int GzipLevel = -1);
    ~DataItemStore();
    DataItemStore(const DataItemStore&) = delete;
    DataItemStore& operator=(const DataItemStore&) = delete;

    /**
     * @brief Store one array and return the text its `<DataItem>` carries.
     * @param rArr The array to store.
     * @return Inline whitespace-separated numbers (`"XML"`), a `.bin` path
     *         (`"Binary"`), or `"<file>.h5:/dataN"` (`"HDF"`).
     * @throws WriteError on an unwritable sibling file, or for `"HDF"` on a
     *         build without HDF5 support.
     */
    std::string Store(const NDArray& rArr);

    /** @brief The data format this store was constructed with. */
    const std::string& DataFormat() const;

    /** @brief Close the HDF5 container, if one was opened. Idempotent. */
    void Close();

    /**
     * @brief Push everything stored so far to disk, keeping the store usable.
     *
     * `"Binary"` writes and closes a file per array, and `"XML"` has no heavy
     * data at all, so this only does real work for `"HDF"` -- but calling it is
     * how a caller says "another process may read this now" without having to
     * know which format it picked.
     */
    void Flush();

    /**
     * @brief The next index the flat `dataN` / `<base>N.bin` naming will use.
     *
     * Exposed only so an appending writer can resume past what an earlier run
     * already wrote; a fresh store starts at 0 and no other caller should care.
     */
    int Counter() const;
    /** @brief Resume the flat naming counter at @p Value (see `Counter()`). */
    void SetCounter(int Value);

    /**
     * @brief Reopen an existing heavy-data container instead of creating one.
     *
     * `"HDF"` opens `<base>.h5` read-write rather than truncating it, which is
     * the difference between continuing a series and destroying it. A missing
     * file is not an error -- appending to a path that does not exist yet is
     * simply a fresh series, which is what makes an append-mode writer usable
     * unconditionally in a restartable solver.
     */
    void OpenExisting();

private:
    struct Impl;
    std::unique_ptr<Impl> mImpl;
};

}  // namespace xdmfcommon
}  // namespace meshioplusplus
