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
// System includes
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/formats/xdmf_time_series.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"

namespace meshioplusplus {

namespace {

// Anonymous-namespace helpers here are `xts_`-prefixed because the single-header
// amalgamation concatenates every src/cpp/src/**.cpp into one translation unit,
// where an unprefixed `add_data_item` would collide with xdmf.cpp's.

/** @brief The XInclude namespace URI the step grids' `xi:include` lives in. */
const char* const xts_xinclude_ns = "http://www.w3.org/2003/XInclude";

/** @brief The name of the single static grid every step's xpointer resolves to. */
const char* const xts_mesh_name = "mesh";

/**
 * @brief Render a time value so it reads back exactly.
 *
 * Shortest of `%.15g`/`%.16g`/`%.17g` that survives a `strtod` round-trip -- an
 * IEEE double always survives at 17, and most step times (`k*dt`) print
 * pleasantly at 15. Deliberately not `std::to_chars`, whose floating-point
 * overload is the documented Emscripten/libc++ hazard this codebase avoids.
 */
std::string xts_format_time(double Value) {
    char buf[40];
    for (int prec = 15; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, Value);
        if (std::strtod(buf, nullptr) == Value)
            return buf;
    }
    return buf;
}

/** @brief Append a `<DataItem>` under @p parent carrying @p rArr via @p rStore. */
void xts_add_data_item(pugi::xml_node parent, xdmfcommon::DataItemStore& rStore,
                       const NDArray& rArr) {
    auto [dtype_s, prec] = xdmfcommon::numpy_to_xdmf_dtype(rArr.Dtype());
    const std::string dims = xdmfcommon::dims_string(rArr);
    pugi::xml_node di = parent.append_child("DataItem");
    di.append_attribute("DataType") = dtype_s;
    di.append_attribute("Dimensions") = dims.c_str();
    di.append_attribute("Format") = rStore.DataFormat().c_str();
    di.append_attribute("Precision") = prec;
    di.text().set(rStore.Store(rArr).c_str());
}

/** @brief Append one `<Attribute>` element (the shape both centers share). */
void xts_add_attribute(pugi::xml_node parent, xdmfcommon::DataItemStore& rStore,
                       const std::string& rName, const char* pCenter, const NDArray& rArr) {
    pugi::xml_node att = parent.append_child("Attribute");
    att.append_attribute("Name") = rName.c_str();
    att.append_attribute("AttributeType") = xdmfcommon::attribute_type(rArr.Shape()).c_str();
    att.append_attribute("Center") = pCenter;
    xts_add_data_item(att, rStore, rArr);
}

}  // namespace

struct XdmfTimeSeriesWriter::Impl {
    std::string mPath;
    pugi::xml_document mDoc;
    pugi::xml_node mCollection;
    std::unique_ptr<xdmfcommon::DataItemStore> mStore;
    std::size_t mNumSteps = 0;
    std::size_t mNumPoints = 0;  // from WritePointsCells, for array validation
    std::size_t mNumCells = 0;
    bool mHasMesh = false;
    bool mFinalized = false;
    bool mAutoFlush = false;

    /**
     * @brief Serialize the document over `mPath` atomically.
     *
     * Written to a sibling temp file and renamed, so a crash mid-write cannot
     * leave a truncated `.xdmf` where a complete earlier one used to be -- the
     * whole point of flushing being that the file on disk is always openable.
     */
    void SaveDocument() {
        const std::string tmp = mPath + ".tmp";
        if (!mDoc.save_file(tmp.c_str(), "  "))
            throw WriteError("XDMF: could not write " + tmp);
        std::error_code ec;
        std::filesystem::rename(tmp, mPath, ec);
        if (ec) {
            std::filesystem::remove(tmp, ec);
            throw WriteError("XDMF: could not write " + mPath);
        }
    }

    /** @brief Start one step's `<Grid>`: the xpointer to the mesh plus `<Time>`. */
    pugi::xml_node BeginStep(double Time) {
        pugi::xml_node grid = mCollection.append_child("Grid");
        // The step references the static grid rather than repeating it. Readers
        // that do not implement XPointer (this repo's own, and Python's
        // TimeSeriesReader) resolve the collection structurally and skip this.
        pugi::xml_node inc = grid.append_child("xi:include");
        const std::string ptr = std::string("xpointer(//Grid[@Name=\"") + xts_mesh_name +
                                "\"]/*[self::Topology or self::Geometry])";
        inc.append_attribute("xpointer") = ptr.c_str();
        pugi::xml_node time = grid.append_child("Time");
        time.append_attribute("Value") = xts_format_time(Time).c_str();
        return grid;
    }
};

XdmfTimeSeriesWriter::XdmfTimeSeriesWriter(const std::string& rPath, const std::string& rDataFormat,
                                           int GzipLevel, XdmfSeriesMode Mode)
    : mImpl(std::make_unique<Impl>()) {
    if (rDataFormat != "XML" && rDataFormat != "Binary" && rDataFormat != "HDF")
        throw WriteError("XDMF: unknown data format '" + rDataFormat +
                         "' (use 'XML', 'Binary', or 'HDF')");
#ifndef MESHIOPLUSPLUS_HAS_HDF5
    if (rDataFormat == "HDF")
        throw WriteError(
            "XDMF: HDF data format requires an HDF5-enabled build "
            "(-DMESHIOPLUSPLUS_WITH_HDF5=ON)");
#endif

    mImpl->mPath = rPath;
    // Sibling heavy-data files, i.e. <path minus extension>.h5 / <...><n>.bin --
    // the same derivation write_xdmf uses, and deliberately not the Python
    // writer's CWD-relative `filename.stem + ".h5"`.
    std::string base = rPath;
    const std::size_t dot = base.find_last_of('.');
    if (dot != std::string::npos)
        base = base.substr(0, dot);
    mImpl->mStore = std::make_unique<xdmfcommon::DataItemStore>(rDataFormat, base, GzipLevel);

    std::error_code ec;
    if (Mode == XdmfSeriesMode::Append && std::filesystem::exists(rPath, ec)) {
        if (!mImpl->mDoc.load_file(rPath.c_str()))
            throw WriteError("XDMF time series: could not parse '" + rPath + "' to append to it");
        // Resolve structurally, the way read_xdmf's xdmf_resolve does -- an
        // xpointer is never followed, and the collection is the first temporal
        // one under <Domain>.
        pugi::xml_node domain = mImpl->mDoc.child("Xdmf").child("Domain");
        for (pugi::xml_node g = domain.child("Grid"); g; g = g.next_sibling("Grid")) {
            const std::string gt = g.attribute("GridType").as_string();
            const std::string ct = g.attribute("CollectionType").as_string();
            if (gt == "Collection" && ct == "Temporal") {
                mImpl->mCollection = g;
                break;
            }
            if (gt == "Uniform" && std::string(g.attribute("Name").as_string()) == xts_mesh_name)
                mImpl->mHasMesh = true;
        }
        if (!mImpl->mCollection)
            throw WriteError("XDMF time series: '" + rPath +
                             "' carries no temporal collection to append to");
        for (pugi::xml_node g = mImpl->mCollection.child("Grid"); g; g = g.next_sibling("Grid"))
            ++mImpl->mNumSteps;
        // A second pass for the mesh grid: it may sit after the collection.
        if (!mImpl->mHasMesh)
            for (pugi::xml_node g = domain.child("Grid"); g; g = g.next_sibling("Grid"))
                if (std::string(g.attribute("Name").as_string()) == xts_mesh_name)
                    mImpl->mHasMesh = true;
        // Resume the heavy-data naming past whatever the earlier run wrote. A
        // mis-resumed counter would silently overwrite data0 rather than fail.
        mImpl->mStore->OpenExisting();
        if (rDataFormat == "Binary") {
            int next = 0;
            while (std::filesystem::exists(base + std::to_string(next) + ".bin", ec))
                ++next;
            mImpl->mStore->SetCounter(next);
        }
        return;
    }

    pugi::xml_node xdmf = mImpl->mDoc.append_child("Xdmf");
    xdmf.append_attribute("Version") = "3.0";
    xdmf.append_attribute("xmlns:xi") = xts_xinclude_ns;
    pugi::xml_node domain = xdmf.append_child("Domain");
    // The collection is created first so the steps accumulate ahead of the static
    // grid in document order, matching the Python writer's layout.
    mImpl->mCollection = domain.append_child("Grid");
    mImpl->mCollection.append_attribute("Name") = "TimeSeries_meshio";
    mImpl->mCollection.append_attribute("GridType") = "Collection";
    mImpl->mCollection.append_attribute("CollectionType") = "Temporal";
}

XdmfTimeSeriesWriter::~XdmfTimeSeriesWriter() {
    if (!mImpl)
        return;  // moved-from
    try {
        Finalize();
    } catch (const std::exception& e) {
        // A destructor must not throw; an explicit Finalize() is how a caller
        // gets to see this failure.
        log::error("XDMF time series: could not finalize '{}': {}", mImpl->mPath, e.what());
    }
}

XdmfTimeSeriesWriter::XdmfTimeSeriesWriter(XdmfTimeSeriesWriter&&) noexcept = default;

XdmfTimeSeriesWriter& XdmfTimeSeriesWriter::operator=(XdmfTimeSeriesWriter&& rOther) noexcept {
    if (this != &rOther) {
        if (mImpl) {
            try {
                Finalize();
            } catch (const std::exception& e) {
                log::error("XDMF time series: could not finalize '{}': {}", mImpl->mPath, e.what());
            }
        }
        mImpl = std::move(rOther.mImpl);
    }
    return *this;
}

void XdmfTimeSeriesWriter::WritePointsCells(const Mesh& rMesh) {
    if (mImpl->mFinalized)
        throw WriteError("XDMF time series: writer already finalized");
    if (mImpl->mHasMesh)
        throw WriteError("XDMF time series: WritePointsCells called more than once");

    pugi::xml_node grid = mImpl->mDoc.child("Xdmf").child("Domain").append_child("Grid");
    grid.append_attribute("Name") = xts_mesh_name;
    grid.append_attribute("GridType") = "Uniform";

    // Geometry
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = points.Shape().size() >= 2 ? points.Shape()[1] : 3;
    if (pdim > 3)
        throw WriteError("XDMF: can only write points up to dimension 3");
    const char* geo_type = (pdim == 1) ? "X" : (pdim == 2) ? "XY" : "XYZ";
    pugi::xml_node geo = grid.append_child("Geometry");
    geo.append_attribute("GeometryType") = geo_type;
    xts_add_data_item(geo, *mImpl->mStore, points);

    mImpl->mNumPoints = rMesh.NumPoints();
    mImpl->mNumCells = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        mImpl->mNumCells += rMesh.Cells(b).NumCells();

    // Topology (single-type, or one packed Mixed array)
    if (rMesh.NumCellBlocks() == 1) {
        const auto cb = rMesh.Cells(0);
        const NDArray& conn = cb.Conn();
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = xdmfcommon::meshio_to_xdmf(cb.Type());
        topo.append_attribute("NumberOfElements") = std::to_string(cb.NumCells()).c_str();
        topo.append_attribute("NodesPerElement") = std::to_string(detail::cols(conn)).c_str();
        xts_add_data_item(topo, *mImpl->mStore, conn);
    } else if (rMesh.NumCellBlocks() > 1) {
        std::size_t total_cells = 0;
        NDArray cd = xdmfcommon::pack_mixed_topology(rMesh, total_cells);
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = "Mixed";
        topo.append_attribute("NumberOfElements") = std::to_string(total_cells).c_str();
        xts_add_data_item(topo, *mImpl->mStore, cd);
    }

    mImpl->mHasMesh = true;
}

void XdmfTimeSeriesWriter::WriteData(double Time, const Mesh& rMesh) {
    if (mImpl->mFinalized)
        throw WriteError("XDMF time series: writer already finalized");
    if (!mImpl->mHasMesh)
        throw WriteError("XDMF time series: WritePointsCells must be called before WriteData");

    pugi::xml_node grid = mImpl->BeginStep(Time);

    // Sorted name order on both, matching write_xdmf and the uniform API's
    // guarantee, so a series is byte-deterministic across mesh backends.
    for (const auto& name : rMesh.PointDataNames())
        xts_add_attribute(grid, *mImpl->mStore, name, "Node", rMesh.PointData(name));

    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        // XDMF stores cell data as one flat array per name across every block.
        NDArray raw = xdmfcommon::concat_cell_data(rMesh, name);
        xts_add_attribute(grid, *mImpl->mStore, name, "Cell", raw);
    }

    ++mImpl->mNumSteps;
    if (mImpl->mAutoFlush)
        Flush();
}

void XdmfTimeSeriesWriter::WriteData(double Time, const std::vector<NamedArray>& rPointData,
                                     const std::vector<NamedArray>& rCellData) {
    if (mImpl->mFinalized)
        throw WriteError("XDMF time series: writer already finalized");
    if (!mImpl->mHasMesh)
        throw WriteError("XDMF time series: WritePointsCells must be called before WriteData");

    pugi::xml_node grid = mImpl->BeginStep(Time);

    const auto emit = [&](const std::vector<NamedArray>& rArrays, std::size_t entities,
                          const char* p_center, const char* p_what) {
        for (const NamedArray& r_a : rArrays) {
            const std::size_t nc = r_a.mNumComponents ? r_a.mNumComponents : 1;
            const std::size_t want = entities * nc;
            if (r_a.mValues.size() != want)
                throw WriteError("XDMF time series: " + std::string(p_what) + " array '" +
                                 r_a.mName + "' holds " + std::to_string(r_a.mValues.size()) +
                                 " values, expected " + std::to_string(want) + " (" +
                                 std::to_string(entities) + " x " + std::to_string(nc) + ")");
            std::vector<std::size_t> shape;
            if (nc == 1)
                shape = {entities};
            else
                shape = {entities, nc};
            NDArray arr = NDArray::Uninit(DType::Float64, shape);
            if (!r_a.mValues.empty())
                std::memcpy(arr.Data(), r_a.mValues.data(), r_a.mValues.size() * sizeof(double));
            // Through the same helper as the Mesh path, so the attribute type,
            // dtype mapping and DataItem shape cannot drift between the two.
            xts_add_attribute(grid, *mImpl->mStore, r_a.mName, p_center, arr);
        }
    };
    emit(rPointData, mImpl->mNumPoints, "Node", "point");
    emit(rCellData, mImpl->mNumCells, "Cell", "cell");

    ++mImpl->mNumSteps;
    if (mImpl->mAutoFlush)
        Flush();
}

void XdmfTimeSeriesWriter::Flush() {
    if (mImpl->mFinalized)
        return;
    // Heavy data first: the .xdmf must never name a dataset that is not on disk.
    mImpl->mStore->Flush();
    mImpl->SaveDocument();
}

void XdmfTimeSeriesWriter::SetAutoFlush(bool Enable) {
    mImpl->mAutoFlush = Enable;
}

bool XdmfTimeSeriesWriter::AutoFlush() const {
    return mImpl && mImpl->mAutoFlush;
}

void XdmfTimeSeriesWriter::Finalize() {
    if (mImpl->mFinalized)
        return;
    // Set first: a failed save must not be retried by the destructor, which
    // could not report the second failure either.
    mImpl->mFinalized = true;
    // Close the heavy data before the light data appears, so a reader that finds
    // the .xdmf always finds a complete .h5 beside it.
    mImpl->mStore->Close();
    mImpl->SaveDocument();
}

std::size_t XdmfTimeSeriesWriter::NumSteps() const {
    return mImpl ? mImpl->mNumSteps : 0;
}

bool XdmfTimeSeriesWriter::Finalized() const {
    return mImpl ? mImpl->mFinalized : true;
}

}  // namespace meshioplusplus
