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

// Compiled only with ADIOS2 (MESHIOPLUSPLUS_HAS_ADIOS2): the whole translation
// unit is empty otherwise, and the registry reports the format as compiled out.
#ifdef MESHIOPLUSPLUS_HAS_ADIOS2

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes
#include <adios2.h>

// Project includes
#include "../detail/library_lock.hpp"
#include "meshioplusplus/formats/vtx.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/vtk_common.hpp"
#include "../detail/typed_view.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "pugixml.hpp"

namespace meshioplusplus {

namespace {

[[noreturn]] void vtx_fail(const std::string& rMessage) {
    throw ReadError("VTX (.bp): " + rMessage);
}

// One data array the schema names: its name in the mesh and the ADIOS2
// variable holding it (the element text when present, else the name).
struct VtxArray {
    std::string mName;
    std::string mVariable;
    bool mCell = false;
};

struct VtxSchema {
    std::string mGeometry = "geometry";
    std::string mConnectivity = "connectivity";
    std::string mTypes = "types";
    std::string mTime;  // empty: no TIME array
    std::vector<VtxArray> mArrays;
};

std::string vtx_trim(std::string S) {
    const char* kWs = " \t\r\n";
    const std::size_t b = S.find_first_not_of(kWs);
    if (b == std::string::npos)
        return {};
    const std::size_t e = S.find_last_not_of(kWs);
    return S.substr(b, e - b + 1);
}

VtxSchema vtx_parse_schema(const std::string& rXml) {
    pugi::xml_document doc;
    if (!doc.load_string(rXml.c_str()))
        vtx_fail("the vtk.xml schema is not valid XML");
    const pugi::xml_node root = doc.child("VTKFile");
    const std::string type = root.attribute("type").as_string();
    if (type != "UnstructuredGrid")
        vtx_fail("only an UnstructuredGrid schema is read, this one is '" + type + "'");
    const pugi::xml_node piece = root.child("UnstructuredGrid").child("Piece");
    if (!piece)
        vtx_fail("the vtk.xml schema has no Piece");
    VtxSchema schema;
    auto variable = [](const pugi::xml_node& rDa) {
        std::string text = vtx_trim(rDa.text().as_string());
        return text.empty() ? std::string(rDa.attribute("Name").as_string()) : text;
    };
    if (const pugi::xml_node da = piece.child("Points").child("DataArray"))
        schema.mGeometry = variable(da);
    for (const pugi::xml_node& rDa : piece.child("Cells").children("DataArray")) {
        const std::string name = rDa.attribute("Name").as_string();
        if (name == "connectivity")
            schema.mConnectivity = variable(rDa);
        else if (name == "types")
            schema.mTypes = variable(rDa);
    }
    for (const char* pGroup : {"PointData", "CellData"}) {
        for (const pugi::xml_node& rDa : piece.child(pGroup).children("DataArray")) {
            const std::string name = rDa.attribute("Name").as_string();
            if (name.empty())
                continue;
            if (name == "TIME") {
                schema.mTime = variable(rDa);
                continue;
            }
            schema.mArrays.push_back({name, variable(rDa), std::string(pGroup) == "CellData"});
        }
    }
    return schema;
}

// What a streaming pass over the steps finds: the schema, the variables each
// step holds, and the step times.
struct VtxLayout {
    VtxSchema mSchema;
    std::vector<std::set<std::string>> mPresent;
    std::vector<double> mTimes;
    bool mShared = false;  // the blocks hold the points they use (see kVtxSharedAttribute)
};

// Set by this library's writer on a file whose blocks each hold the points they use (several
// cell types, or the parts of a partition): the ghost copies are one mesh and are welded on read.
constexpr const char* kVtxSharedAttribute = "meshio.shared_points";

// A single value of the current step. A value several ranks wrote has one
// instance per writer, so it is taken from the block list: a Get would deliver
// every instance into a one-element buffer.
template <class T>
double vtx_first_value(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName) {
    adios2::Variable<T> var = rIo.InquireVariable<T>(rName);
    const auto infos = rEngine.BlocksInfo(var, rEngine.CurrentStep());
    if (infos.empty())
        vtx_fail("the time variable '" + rName + "' holds no value");
    if (infos.front().IsValue)
        return static_cast<double>(infos.front().Value);
    std::vector<T> v(1);
    var.SetBlockSelection(0);
    rEngine.Get(var, v, adios2::Mode::Sync);
    return v.empty() ? 0.0 : static_cast<double>(v.front());
}

double vtx_get_scalar(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName) {
    const std::string type = rIo.VariableType(rName);
    if (type == "double")
        return vtx_first_value<double>(rIo, rEngine, rName);
    if (type == "float")
        return vtx_first_value<float>(rIo, rEngine, rName);
    vtx_fail("the time variable '" + rName + "' is of type " + type);
}

VtxLayout vtx_scan(const std::string& rPath) {
    VtxLayout layout;
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("meshioplusplus_vtx_scan");
    adios2::Engine engine = io.Open(rPath, adios2::Mode::Read);
    bool have_schema = false;
    while (engine.BeginStep() == adios2::StepStatus::OK) {
        if (!have_schema) {
            const adios2::Attribute<std::string> attr = io.InquireAttribute<std::string>("vtk.xml");
            if (!attr) {
                engine.Close();
                vtx_fail("'" + rPath +
                         "' has no vtk.xml attribute: only DOLFINx VTXWriter output is read "
                         "(adios4dolfinx checkpoints and Fides output are other layouts)");
            }
            std::string xml;
            for (const std::string& rPart : attr.Data())
                xml += rPart;
            layout.mSchema = vtx_parse_schema(xml);
            layout.mShared =
                static_cast<bool>(io.InquireAttribute<std::string>(kVtxSharedAttribute));
            have_schema = true;
        }
        std::set<std::string> names;
        for (const auto& rEntry : io.AvailableVariables())
            names.insert(rEntry.first);
        const std::size_t k = layout.mPresent.size();
        const std::string& rTime = layout.mSchema.mTime;
        layout.mTimes.push_back(!rTime.empty() && names.count(rTime)
                                    ? vtx_get_scalar(io, engine, rTime)
                                    : static_cast<double>(k));
        layout.mPresent.push_back(std::move(names));
        engine.EndStep();
    }
    engine.Close();
    if (!have_schema)
        vtx_fail("'" + rPath + "' holds no step");
    return layout;
}

// A variable's index among the steps that hold it: ADIOS2's random-access
// step selection and block lists count only those.
std::size_t vtx_relative_step(const VtxLayout& rLayout, const std::string& rName,
                              std::size_t Step) {
    std::size_t r = 0;
    for (std::size_t s = 0; s < Step; ++s)
        r += rLayout.mPresent[s].count(rName);
    return r;
}

// The last step at or before `Step` holding `rName`, or npos.
std::size_t vtx_last_with(const VtxLayout& rLayout, const std::string& rName, std::size_t Step) {
    for (std::size_t s = Step + 1; s-- > 0;)
        if (rLayout.mPresent[s].count(rName))
            return s;
    return static_cast<std::size_t>(-1);
}

template <class T>
DType vtx_dtype();
template <>
DType vtx_dtype<double>() {
    return DType::Float64;
}
template <>
DType vtx_dtype<float>() {
    return DType::Float32;
}
template <>
DType vtx_dtype<std::int8_t>() {
    return DType::Int8;
}
template <>
DType vtx_dtype<std::int16_t>() {
    return DType::Int16;
}
template <>
DType vtx_dtype<std::int32_t>() {
    return DType::Int32;
}
template <>
DType vtx_dtype<std::int64_t>() {
    return DType::Int64;
}
template <>
DType vtx_dtype<std::uint8_t>() {
    return DType::UInt8;
}
template <>
DType vtx_dtype<std::uint16_t>() {
    return DType::UInt16;
}
template <>
DType vtx_dtype<std::uint32_t>() {
    return DType::UInt32;
}
template <>
DType vtx_dtype<std::uint64_t>() {
    return DType::UInt64;
}

// Every block of a variable at one (relative) step, each as an NDArray shaped
// by the block's count (a value block is one element). With `pRows`, only the
// blocks matched to the rank blocks by length, in order, are read: DOLFINx
// 0.11 writes stray blocks whose count does not describe their payload (ADIOS2
// then reads past the buffer it sized), so an unmatched block is never read.
template <class T>
std::vector<NDArray> vtx_blocks_typed(adios2::IO& rIo, adios2::Engine& rEngine,
                                      const std::string& rName, std::size_t Rel,
                                      const std::vector<std::size_t>* pRows) {
    adios2::Variable<T> var = rIo.InquireVariable<T>(rName);
    std::vector<NDArray> out;
    if (!var)
        return out;
    const auto infos = rEngine.BlocksInfo(var, Rel);
    var.SetStepSelection({Rel, 1});
    if (infos.empty() && var.ShapeID() == adios2::ShapeID::GlobalValue) {
        std::vector<T> values;
        rEngine.Get(var, values, adios2::Mode::Sync);
        if (!values.empty()) {
            NDArray a(vtx_dtype<T>(), {1});
            a.As<T>()[0] = values.front();
            out.push_back(std::move(a));
        }
        return out;
    }
    std::vector<bool> wanted(infos.size(), pRows == nullptr);
    if (pRows) {
        std::size_t j = 0;
        for (std::size_t b = 0; b < infos.size(); ++b) {
            const std::size_t rows =
                infos[b].IsValue ? 1 : (infos[b].Count.empty() ? 1 : infos[b].Count[0]);
            if (j < pRows->size() && rows == (*pRows)[j]) {
                wanted[b] = true;
                ++j;
            } else {
                log::debug("VTX: skipping block {} of '{}': it matches no rank's size", b, rName);
            }
        }
    }
    // A global value several ranks wrote lists one instance per writer, each
    // carrying its value: every block below is then an IsValue block.
    for (std::size_t b = 0; b < infos.size(); ++b) {
        if (!wanted[b])
            continue;
        std::vector<std::size_t> shape;
        if (infos[b].IsValue) {
            NDArray a(vtx_dtype<T>(), {1});
            a.As<T>()[0] = infos[b].Value;
            out.push_back(std::move(a));
            continue;
        }
        for (const std::size_t c : infos[b].Count)
            shape.push_back(c);
        if (shape.empty())
            shape.push_back(1);
        NDArray a(vtx_dtype<T>(), shape);
        if (a.Size() > 0) {
            // Into a vector ADIOS2 sizes itself, then checked: a raw buffer
            // would be overrun by any selection larger than the block list says.
            std::vector<T> values;
            var.SetBlockSelection(b);
            rEngine.Get(var, values, adios2::Mode::Sync);
            if (values.size() != a.Size())
                vtx_fail("block " + std::to_string(b) + " of '" + rName + "' holds " +
                         std::to_string(values.size()) + " values, its count " +
                         std::to_string(a.Size()));
            std::copy(values.begin(), values.end(), a.As<T>());
        }
        out.push_back(std::move(a));
    }
    return out;
}

std::vector<NDArray> vtx_blocks(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName,
                                std::size_t Rel, const std::vector<std::size_t>* pRows = nullptr) {
    const std::string type = rIo.VariableType(rName);
    if (type == "double")
        return vtx_blocks_typed<double>(rIo, rEngine, rName, Rel, pRows);
    if (type == "float")
        return vtx_blocks_typed<float>(rIo, rEngine, rName, Rel, pRows);
    if (type == "int8_t")
        return vtx_blocks_typed<std::int8_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "int16_t")
        return vtx_blocks_typed<std::int16_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "int32_t")
        return vtx_blocks_typed<std::int32_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "int64_t")
        return vtx_blocks_typed<std::int64_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "uint8_t")
        return vtx_blocks_typed<std::uint8_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "uint16_t")
        return vtx_blocks_typed<std::uint16_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "uint32_t")
        return vtx_blocks_typed<std::uint32_t>(rIo, rEngine, rName, Rel, pRows);
    if (type == "uint64_t")
        return vtx_blocks_typed<std::uint64_t>(rIo, rEngine, rName, Rel, pRows);
    if (type.empty())
        return {};
    vtx_fail("variable '" + rName + "' has the unsupported type " + type);
}

// The row count of every block of a variable at one (relative) step, from the
// block list alone (a value block counts one row).
std::vector<std::size_t> vtx_block_rows(adios2::IO& rIo, adios2::Engine& rEngine,
                                        const std::string& rName, std::size_t Rel) {
    std::vector<std::size_t> rows;
    auto collect = [&](auto Tag) {
        using T = decltype(Tag);
        adios2::Variable<T> var = rIo.InquireVariable<T>(rName);
        if (!var)
            return;
        for (const auto& rInfo : rEngine.BlocksInfo(var, Rel))
            rows.push_back(rInfo.IsValue || rInfo.Count.empty() ? 1 : rInfo.Count[0]);
    };
    const std::string type = rIo.VariableType(rName);
    if (type == "double")
        collect(double{});
    else if (type == "float")
        collect(float{});
    else if (type == "int8_t")
        collect(std::int8_t{});
    else if (type == "int16_t")
        collect(std::int16_t{});
    else if (type == "int32_t")
        collect(std::int32_t{});
    else if (type == "int64_t")
        collect(std::int64_t{});
    else if (type == "uint8_t")
        collect(std::uint8_t{});
    else if (type == "uint16_t")
        collect(std::uint16_t{});
    else if (type == "uint32_t")
        collect(std::uint32_t{});
    else if (type == "uint64_t")
        collect(std::uint64_t{});
    return rows;
}

std::int64_t vtx_int_at(const NDArray& rA, std::size_t i) {
    switch (rA.Dtype()) {
        case DType::Int8:
            return rA.As<std::int8_t>()[i];
        case DType::Int16:
            return rA.As<std::int16_t>()[i];
        case DType::Int32:
            return rA.As<std::int32_t>()[i];
        case DType::Int64:
            return rA.As<std::int64_t>()[i];
        case DType::UInt8:
            return rA.As<std::uint8_t>()[i];
        case DType::UInt16:
            return rA.As<std::uint16_t>()[i];
        case DType::UInt32:
            return rA.As<std::uint32_t>()[i];
        case DType::UInt64:
            return static_cast<std::int64_t>(rA.As<std::uint64_t>()[i]);
        default:
            vtx_fail("an integer variable holds floating-point values");
    }
}

// A VTK_LAGRANGE_* cell whose node order is a linear or quadratic VTK cell's
// is read as that cell (see vtx.hpp); anything else keeps its type.
std::int64_t vtx_lower_type(std::int64_t Type, std::size_t Nodes) {
    struct Rule {
        std::int64_t mLagrange;
        std::size_t mNodes;
        std::int64_t mVtk;
    };
    static const Rule kRules[] = {
        {68, 2, 3},  {68, 3, 21},  {69, 3, 5},  {69, 6, 22}, {70, 4, 9},  {70, 9, 28},
        {71, 4, 10}, {71, 10, 24}, {72, 8, 12}, {73, 6, 13}, {74, 5, 14},
    };
    for (const Rule& rR : kRules)
        if (rR.mLagrange == Type && rR.mNodes == Nodes)
            return rR.mVtk;
    return Type;
}

std::size_t vtx_rows(const NDArray& rA) {
    return rA.Shape().empty() ? 1 : rA.Shape()[0];
}

std::size_t vtx_row_width(const NDArray& rA) {
    std::size_t w = 1;
    for (std::size_t d = 1; d < rA.Shape().size(); ++d)
        w *= rA.Shape()[d];
    return w;
}

// Match an array's blocks to the rank blocks by length, in order; a block of
// any other length is skipped. Returns the concatenation, or an empty array
// when the blocks do not cover every rank.
NDArray vtx_concat_matching(const std::vector<NDArray>& rBlocks,
                            const std::vector<std::size_t>& rCounts, const std::string& rName) {
    std::vector<const NDArray*> picked;
    std::size_t j = 0;
    for (const NDArray& rB : rBlocks) {
        if (j < rCounts.size() && vtx_rows(rB) == rCounts[j]) {
            picked.push_back(&rB);
            ++j;
        } else {
            log::debug("VTX: skipping a block of '{}' that matches no rank's size", rName);
        }
    }
    if (j != rCounts.size() || picked.empty())
        return {};
    const std::size_t width = vtx_row_width(*picked.front());
    std::size_t total = 0;
    for (const NDArray* pB : picked) {
        if (vtx_row_width(*pB) != width || pB->Dtype() != picked.front()->Dtype())
            return {};
        total += vtx_rows(*pB);
    }
    std::vector<std::size_t> shape = {total};
    if (picked.front()->Shape().size() > 1 && width > 1)
        shape.push_back(width);
    NDArray out = NDArray::Uninit(picked.front()->Dtype(), shape);
    const std::size_t item = dtype_size(out.Dtype()) * width;
    std::byte* pDst = out.Data();
    for (const NDArray* pB : picked) {
        const std::size_t n = vtx_rows(*pB) * item;
        if (n > 0)
            std::memcpy(pDst, pB->Data(), n);
        pDst += n;
    }
    return out;
}

// Keep the rows `rKeep` of a per-point array (ghost welding).
NDArray vtx_take_rows(const NDArray& rA, const std::vector<std::size_t>& rKeep) {
    std::vector<std::size_t> shape = rA.Shape();
    shape[0] = rKeep.size();
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);
    const std::size_t item = dtype_size(rA.Dtype()) * vtx_row_width(rA);
    for (std::size_t i = 0; i < rKeep.size(); ++i)
        std::memcpy(out.Data() + i * item, rA.Data() + rKeep[i] * item, item);
    return out;
}

struct VtxStep {
    NDArray mPoints;
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mOffsets;
    std::vector<std::int64_t> mTypes;
    std::map<std::string, NDArray> mPointData;
    std::unordered_map<std::string, NDArray> mCellData;
    double mTime = 0.0;
};

// GhostPolicy::Drop: weld every point onto the owned copy sharing its
// vtkOriginalPointIds (the first with vtkGhostType 0, else the first seen), keep
// the surviving points in their first-seen order, and drop both ghost arrays.
void vtx_weld_ghosts(VtxStep& rStep, bool SortById) {
    auto ids_it = rStep.mPointData.find("vtkOriginalPointIds");
    auto ghost_it = rStep.mPointData.find("vtkGhostType");
    if (ids_it == rStep.mPointData.end()) {
        log::warn("VTX: no vtkOriginalPointIds, so ghost points cannot be welded");
        return;
    }
    const NDArray& rIds = ids_it->second;
    const std::size_t n = vtx_rows(rIds);
    std::unordered_map<std::int64_t, std::size_t> owner;
    owner.reserve(n);
    for (int pass = 0; pass < 2; ++pass)
        for (std::size_t i = 0; i < n; ++i) {
            const bool ghost =
                ghost_it != rStep.mPointData.end() && vtx_int_at(ghost_it->second, i) != 0;
            if ((pass == 0) != ghost)
                owner.emplace(vtx_int_at(rIds, i), i);
        }
    std::vector<std::size_t> keep;
    std::vector<std::int64_t> new_index(n, -1);
    for (std::size_t i = 0; i < n; ++i)
        if (owner[vtx_int_at(rIds, i)] == i)
            keep.push_back(i);
    // A shared-points file lists the points in the mesh's order, whatever block holds them.
    if (SortById)
        std::stable_sort(keep.begin(), keep.end(), [&](std::size_t A, std::size_t B) {
            return vtx_int_at(rIds, A) < vtx_int_at(rIds, B);
        });
    for (std::size_t k = 0; k < keep.size(); ++k)
        new_index[keep[k]] = static_cast<std::int64_t>(k);
    for (std::int64_t& rV : rStep.mConn)
        rV = new_index[owner[vtx_int_at(rIds, static_cast<std::size_t>(rV))]];
    rStep.mPoints = vtx_take_rows(rStep.mPoints, keep);
    rStep.mPointData.erase("vtkGhostType");
    rStep.mPointData.erase("vtkOriginalPointIds");
    for (auto& rEntry : rStep.mPointData)
        rEntry.second = vtx_take_rows(rEntry.second, keep);
}

VtxStep vtx_read_step(const std::string& rPath, const ReadOptions& rOpts) {
    const VtxLayout layout = vtx_scan(rPath);
    const VtxSchema& rS = layout.mSchema;
    const std::size_t step = rOpts.ResolveTimeStep(layout.mPresent.size());
    const std::size_t mesh_step = vtx_last_with(layout, rS.mGeometry, step);
    if (mesh_step == static_cast<std::size_t>(-1))
        vtx_fail("no step up to " + std::to_string(step) + " holds the mesh ('" + rS.mGeometry +
                 "')");
    const std::size_t conn_step = vtx_last_with(layout, rS.mConnectivity, step);
    const std::size_t types_step = vtx_last_with(layout, rS.mTypes, step);
    if (conn_step == static_cast<std::size_t>(-1) || types_step == static_cast<std::size_t>(-1))
        vtx_fail("no step up to " + std::to_string(step) + " holds the cells");

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("meshioplusplus_vtx_read");
    adios2::Engine engine = io.Open(rPath, adios2::Mode::ReadRandomAccess);
    auto blocks_at = [&](const std::string& rName, std::size_t S,
                         const std::vector<std::size_t>* pRows = nullptr) {
        return vtx_blocks(io, engine, rName, vtx_relative_step(layout, rName, S), pRows);
    };

    VtxStep out;
    out.mTime = layout.mTimes[step];

    // Points: one block per rank.
    const std::vector<NDArray> geom = blocks_at(rS.mGeometry, mesh_step);
    std::vector<std::size_t> point_counts;
    for (const NDArray& rG : geom)
        point_counts.push_back(vtx_rows(rG));
    out.mPoints = vtx_concat_matching(geom, point_counts, rS.mGeometry);
    if (out.mPoints.Shape().size() != 2)
        vtx_fail("'" + rS.mGeometry + "' is not an (n, dim) array");
    if (out.mPoints.Dtype() != DType::Float64 && out.mPoints.Dtype() != DType::Float32)
        vtx_fail("'" + rS.mGeometry + "' is not floating point");

    // Cells: rank blocks of (cells, 1 + nodes) rows, local point numbering.
    const std::vector<NDArray> conn = blocks_at(rS.mConnectivity, conn_step);
    if (conn.size() != geom.size())
        vtx_fail("'" + rS.mConnectivity + "' has " + std::to_string(conn.size()) + " blocks, '" +
                 rS.mGeometry + "' " + std::to_string(geom.size()));
    const std::vector<NDArray> types = blocks_at(rS.mTypes, types_step);
    if (types.empty())
        vtx_fail("'" + rS.mTypes + "' holds no value");
    std::vector<std::size_t> cell_counts;
    std::int64_t base = 0;
    std::int64_t end = 0;
    for (std::size_t b = 0; b < conn.size(); ++b) {
        const NDArray& rC = conn[b];
        const std::size_t rows = rC.Shape().size() == 2 ? rC.Shape()[0] : 0;
        const std::size_t width = rC.Shape().size() == 2 ? rC.Shape()[1] : 0;
        if (rows > 0 && width < 2)
            vtx_fail("'" + rS.mConnectivity + "' is not a (cells, 1 + nodes) array");
        cell_counts.push_back(rows);
        // `types` is one global value, one value per block, or one per cell.
        const NDArray& rT = types.size() == conn.size() ? types[b] : types.front();
        for (std::size_t c = 0; c < rows; ++c) {
            const std::int64_t n = vtx_int_at(rC, c * width);
            if (n != static_cast<std::int64_t>(width) - 1)
                vtx_fail("a row of '" + rS.mConnectivity + "' counts " + std::to_string(n) +
                         " nodes in a row of " + std::to_string(width - 1));
            for (std::size_t j = 1; j < width; ++j) {
                const std::int64_t v = vtx_int_at(rC, c * width + j);
                if (v < 0 || static_cast<std::size_t>(v) >= point_counts[b])
                    vtx_fail("'" + rS.mConnectivity + "' names point " + std::to_string(v) +
                             " of a block with " + std::to_string(point_counts[b]));
                out.mConn.push_back(v + base);
            }
            end += static_cast<std::int64_t>(width) - 1;
            out.mOffsets.push_back(end);
            const std::size_t ti = rT.Size() == rows ? c : 0;
            out.mTypes.push_back(vtx_lower_type(vtx_int_at(rT, ti), width - 1));
        }
        base += static_cast<std::int64_t>(point_counts[b]);
    }

    // Data arrays of the requested step itself.
    const bool drop = rOpts.mGhosts == GhostPolicy::Drop || layout.mShared;
    for (const VtxArray& rA : rS.mArrays) {
        const bool ghost_array = rA.mName == "vtkGhostType" || rA.mName == "vtkOriginalPointIds";
        const bool wanted = !rOpts.mPointsOnly && rOpts.WantsArray(rA.mName);
        if (!wanted && !(drop && ghost_array && !rA.mCell))
            continue;
        // The ghost arrays of a `reuse` file are written with the mesh only.
        constexpr std::size_t kNone = static_cast<std::size_t>(-1);
        std::size_t start = step;
        if (!layout.mPresent[step].count(rA.mVariable)) {
            start = ghost_array ? vtx_last_with(layout, rA.mVariable, step) : kNone;
            if (start == kNone)
                continue;
        }
        const std::vector<std::size_t>& rRows = rA.mCell ? cell_counts : point_counts;
        // DOLFINx 0.11 writes the ghost arrays twice in the first step of a
        // `reuse` file, the second time as one-value blocks whose payload is not
        // what their count says -- and ADIOS2 then misreads the next rank's
        // block too. An array whose blocks do not match the rank blocks one to
        // one is taken from the nearest step where they do (DOLFINx rewrites
        // the ghost arrays, unchanged, every step).
        std::size_t source = start;
        auto clean = [&](std::size_t S) {
            return layout.mPresent[S].count(rA.mVariable) &&
                   vtx_block_rows(io, engine, rA.mVariable,
                                  vtx_relative_step(layout, rA.mVariable, S)) == rRows;
        };
        if (!clean(start)) {
            source = kNone;
            // Only the ghost arrays are constant over the steps; any other array
            // is dropped rather than replaced by another step's values.
            for (std::size_t d = 1; ghost_array && d < layout.mPresent.size() && source == kNone;
                 ++d) {
                if (start + d < layout.mPresent.size() && clean(start + d))
                    source = start + d;
                else if (d <= start && clean(start - d))
                    source = start - d;
            }
            if (source == kNone) {
                log::warn("VTX: the blocks of '{}' match no rank's {} count and it is not read",
                          rA.mName, rA.mCell ? "cell" : "point");
                continue;
            }
            log::debug("VTX: '{}' taken from step {}: its blocks at step {} are malformed",
                       rA.mName, source, step);
        }
        const std::vector<NDArray> blocks = blocks_at(rA.mVariable, source, &rRows);
        NDArray joined = vtx_concat_matching(blocks, rRows, rA.mName);
        if (joined.Shape().empty()) {
            log::warn("VTX: '{}' matches no rank's {} count and is not read", rA.mName,
                      rA.mCell ? "cell" : "point");
            continue;
        }
        if (rA.mCell)
            out.mCellData.emplace(rA.mName, std::move(joined));
        else
            out.mPointData.emplace(rA.mName, std::move(joined));
    }
    engine.Close();

    if (drop)
        vtx_weld_ghosts(out, layout.mShared);
    return out;
}

}  // namespace

Mesh read_vtx(const std::string& rPath, const ReadOptions& rOpts) {
    detail::LibraryLock lock;
    VtxStep step;
    try {
        step = vtx_read_step(rPath, rOpts);
    } catch (const ReadError&) {
        throw;
    } catch (const std::exception& rE) {
        // ADIOS2 reports a missing or foreign file with its own exceptions.
        vtx_fail(std::string(rE.what()));
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(step.mPoints));
    detail::check_vtk_cell_arrays(step.mConn.size(), step.mOffsets, step.mTypes, step.mCellData);
    detail::reconstruct_cells(step.mConn.data(), step.mOffsets, step.mTypes, step.mCellData, mesh);
    for (auto& rEntry : step.mPointData)
        mesh.AddPointData(rEntry.first, std::move(rEntry.second));
    NDArray t(DType::Float64, {1});
    t.As<double>()[0] = step.mTime;
    mesh.AddFieldData(kSequenceTimeKey, std::move(t));
    return mesh;
}

std::vector<double> vtx_time_values(const std::string& rPath) {
    detail::LibraryLock lock;
    try {
        return vtx_scan(rPath).mTimes;
    } catch (const ReadError&) {
        throw;
    } catch (const std::exception& rE) {
        vtx_fail(std::string(rE.what()));
    }
}

MeshMetadata read_vtx_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    detail::LibraryLock lock;
    if (rOpts.mGhosts == GhostPolicy::Drop) {
        // The welded point count needs the ids: read the step.
        ReadOptions options = rOpts;
        options.mTimeStep = 0;
        MeshMetadata meta = metadata_from_mesh(read_vtx(rPath, options));
        meta.mFellBackToFullRead = true;
        meta.mFormat = "vtx";
        meta.mTimeValues = vtx_time_values(rPath);
        return meta;
    }
    MeshMetadata meta;
    meta.mFormat = "vtx";
    try {
        const VtxLayout layout = vtx_scan(rPath);
        const VtxSchema& rS = layout.mSchema;
        meta.mTimeValues = layout.mTimes;
        if (!layout.mPresent[0].count(rS.mGeometry) || !layout.mPresent[0].count(rS.mTypes) ||
            !layout.mPresent[0].count(rS.mConnectivity))
            vtx_fail("step 0 holds no mesh");
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("meshioplusplus_vtx_meta");
        adios2::Engine engine = io.Open(rPath, adios2::Mode::ReadRandomAccess);
        // Block sizes from the metadata; `types` is a handful of values.
        auto counts = [&](const std::string& rName) {
            std::vector<std::vector<std::size_t>> out;
            const std::string type = io.VariableType(rName);
            auto collect = [&](auto Tag) {
                using T = decltype(Tag);
                adios2::Variable<T> var = io.InquireVariable<T>(rName);
                for (const auto& rInfo : engine.BlocksInfo(var, 0)) {
                    std::vector<std::size_t> c(rInfo.Count.begin(), rInfo.Count.end());
                    out.push_back(c);
                }
            };
            if (type == "double")
                collect(double{});
            else if (type == "float")
                collect(float{});
            else if (type == "int64_t")
                collect(std::int64_t{});
            else if (type == "int32_t")
                collect(std::int32_t{});
            else
                vtx_fail("variable '" + rName + "' has the unsupported type " + type);
            return out;
        };
        std::size_t npoints = 0;
        for (const auto& rC : counts(rS.mGeometry)) {
            npoints += rC.empty() ? 0 : rC[0];
            meta.mPointDim = rC.size() > 1 ? rC[1] : 3;
        }
        meta.mNumPoints = npoints;
        const std::vector<NDArray> types = vtx_blocks(io, engine, rS.mTypes, 0);
        const auto conn = counts(rS.mConnectivity);
        std::vector<std::int64_t> cell_types;
        std::vector<std::int64_t> offsets;
        std::int64_t end = 0;
        for (std::size_t b = 0; b < conn.size(); ++b) {
            const std::size_t rows = conn[b].size() == 2 ? conn[b][0] : 0;
            const std::size_t width = conn[b].size() == 2 ? conn[b][1] : 0;
            const NDArray& rT = types.size() == conn.size() ? types[b] : types.front();
            for (std::size_t c = 0; c < rows; ++c) {
                end += static_cast<std::int64_t>(width) - 1;
                offsets.push_back(end);
                cell_types.push_back(
                    vtx_lower_type(vtx_int_at(rT, rT.Size() == rows ? c : 0), width - 1));
            }
        }
        engine.Close();
        meta.mCellBlocks = detail::summarize_cells(offsets, cell_types);
        if (!rOpts.mPointsOnly)
            for (const VtxArray& rA : rS.mArrays)
                if (layout.mPresent[0].count(rA.mVariable) && rOpts.WantsArray(rA.mName))
                    (rA.mCell ? meta.mCellDataNames : meta.mPointDataNames).push_back(rA.mName);
        std::sort(meta.mPointDataNames.begin(), meta.mPointDataNames.end());
        std::sort(meta.mCellDataNames.begin(), meta.mCellDataNames.end());
        meta.mFieldDataNames.push_back(kSequenceTimeKey);
    } catch (const ReadError&) {
        throw;
    } catch (const std::exception& rE) {
        vtx_fail(std::string(rE.what()));
    }
    return meta;
}

// ---------------------------------------------------------------------------
// The writer: VTXWriter's layout. One block for a mesh of one cell type, as for one rank;
// a mesh of several cell types, or one whose cells carry `partition:part` values of several
// parts, is one block per cell type (and part), every block holding the points it uses. The
// Python twin (`vtx/_write.py`) writes the same variables, attributes and values.
// ---------------------------------------------------------------------------
namespace {

constexpr const char* kVwPartition = "partition:part";
constexpr const char* kVwGhostType = "vtkGhostType";
constexpr const char* kVwOriginalIds = "vtkOriginalPointIds";

[[noreturn]] void vw_fail(const std::string& rMessage) {
    throw WriteError("VTX (.bp): " + rMessage);
}

// (VTK type, node count) -> VTK_LAGRANGE_* type, the inverse of the reader's lowering.
int vw_lagrange(int VtkType, std::size_t Nodes) {
    static const std::map<std::pair<int, std::size_t>, int> kMap = {
        {{3, 2}, 68},  {{21, 3}, 68},  {{5, 3}, 69},  {{22, 6}, 69}, {{9, 4}, 70}, {{28, 9}, 70},
        {{10, 4}, 71}, {{24, 10}, 71}, {{12, 8}, 72}, {{13, 6}, 73}, {{14, 5}, 74}};
    const auto it = kMap.find({VtkType, Nodes});
    return it == kMap.end() ? VtkType : it->second;
}

bool vw_reserved(const std::string& rName) {
    return rName == "geometry" || rName == "connectivity" || rName == "types" || rName == "step" ||
           rName == "TIME" || rName == "NumberOfNodes" || rName == "NumberOfEntities";
}

// Python's xml.sax.saxutils.quoteattr, so the schema attribute reads the same.
std::string vw_quoteattr(const std::string& rValue) {
    std::string v;
    for (const char c : rValue) {
        if (c == '&')
            v += "&amp;";
        else if (c == '<')
            v += "&lt;";
        else if (c == '>')
            v += "&gt;";
        else if (c == '\n')
            v += "&#10;";
        else if (c == '\r')
            v += "&#13;";
        else if (c == '\t')
            v += "&#9;";
        else
            v += c;
    }
    char quote = '"';
    if (v.find('"') != std::string::npos) {
        if (v.find('\'') != std::string::npos) {
            std::string e;
            for (const char c : v)
                e += c == '"' ? std::string("&quot;") : std::string(1, c);
            v = e;
        } else {
            quote = '\'';
        }
    }
    return std::string(1, quote) + v + std::string(1, quote);
}

struct VwArray {
    std::string mName;
    bool mCell = false;
};

// One block of the file: a cell type (of one part) and the points it uses.
struct VwGroup {
    int mVtkType = 0;
    std::size_t mWidth = 0;             // 1 + nodes
    std::vector<std::int64_t> mConn;    // rows x mWidth, each row led by its node count
    bool mAllCells = true;              // the block's cells are the mesh's, in order
    std::vector<std::size_t> mCells;    // else their numbers in the mesh
    bool mAllPoints = true;             // the block's points are the mesh's
    std::vector<std::int64_t> mPoints;  // else the mesh points it holds
    std::vector<std::uint8_t> mGhost;   // 1 where another block owns the point
    std::size_t Rows() const { return mWidth == 0 ? 0 : mConn.size() / mWidth; }
};

// The rows `rIndex` of an array of `Rows` rows, as bytes.
std::vector<std::byte> vw_take_rows(const NDArray& rArray, std::size_t Rows,
                                    const std::vector<std::size_t>& rIndex) {
    const std::size_t row_bytes = Rows == 0 ? 0 : rArray.Size() / Rows * dtype_size(rArray.Dtype());
    const std::byte* p = static_cast<const std::byte*>(rArray.Data());
    std::vector<std::byte> out;
    out.reserve(rIndex.size() * row_bytes);
    for (const std::size_t i : rIndex)
        out.insert(out.end(), p + i * row_bytes, p + (i + 1) * row_bytes);
    return out;
}

class VwModel {
public:
    explicit VwModel(const Mesh& rMesh) {
        mN = rMesh.NumPoints();
        mGeometry.assign(mN * 3, 0.0);
        {
            const NDArray& pts = rMesh.Points();
            const detail::DoubleView view(pts);
            const std::size_t width = mN == 0 ? 3 : pts.Size() / mN;
            for (std::size_t i = 0; i < mN; ++i)
                for (std::size_t k = 0; k < 3 && k < width; ++k)
                    mGeometry[i * 3 + k] = view.Data()[i * width + k];
        }
        struct Block {
            std::size_t mIndex = 0;
            std::string mType;
            int mVtk = 0;
            std::size_t mWidth = 0, mFirst = 0, mCount = 0;
            std::vector<std::int64_t> mRows;
        };
        std::vector<Block> blocks;
        const auto& vtk_types = meshio_to_vtk_type();
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const std::size_t count = cb.NumCells();
            if (count == 0)
                continue;
            const auto found = vtk_types.find(cb.Type());
            if (found == vtk_types.end())
                vw_fail("cell type '" + cb.Type() + "' has no VTK type");
            const std::size_t w = cb.NodesPerCell();
            std::vector<int> order = meshio_to_vtk_order(cb.Type());
            if (order.size() != w)
                order.clear();
            Block blk;
            blk.mIndex = b;
            blk.mType = cb.Type();
            blk.mVtk = vw_lagrange(found->second, w);
            blk.mWidth = w + 1;
            blk.mFirst = mCells;
            blk.mCount = count;
            const detail::Int64View conn(cb.Conn());
            blk.mRows.reserve(count * (w + 1));
            for (std::size_t r = 0; r < count; ++r) {
                blk.mRows.push_back(static_cast<std::int64_t>(w));
                for (std::size_t j = 0; j < w; ++j) {
                    const std::int64_t node =
                        conn.Data()[r * w +
                                    (order.empty() ? j : static_cast<std::size_t>(order[j]))];
                    if (node < 0 || static_cast<std::size_t>(node) >= mN)
                        vw_fail("a '" + cb.Type() +
                                "' cell refers to a point the mesh does not have");
                    blk.mRows.push_back(node);
                }
            }
            mCells += count;
            blocks.push_back(std::move(blk));
        }
        if (mCells == 0)
            vw_fail("the mesh has no cells");

        // The parts, when every cell block says which one each cell belongs to.
        std::vector<std::vector<std::int64_t>> parts;
        if (rMesh.HasCellData(kVwPartition) &&
            rMesh.CellDataNumBlocks(kVwPartition) == rMesh.NumCellBlocks()) {
            bool ok = true;
            for (const Block& blk : blocks) {
                const NDArray& a = rMesh.CellData(kVwPartition, blk.mIndex);
                if (detail::is_float_dtype(a.Dtype()) || a.Shape().size() != 1 ||
                    a.Size() != blk.mCount) {
                    ok = false;
                    break;
                }
                const detail::Int64View view(a);
                parts.emplace_back(view.Data(), view.Data() + blk.mCount);
            }
            if (!ok)
                parts.clear();
        }
        std::set<std::int64_t> labels;
        for (const auto& rPart : parts)
            labels.insert(rPart.begin(), rPart.end());
        std::set<std::string> types;
        for (const Block& blk : blocks)
            types.insert(blk.mType);

        auto from_rows = [](const Block& rBlk, const std::vector<std::size_t>* pKeep) {
            VwGroup g;
            g.mVtkType = rBlk.mVtk;
            g.mWidth = rBlk.mWidth;
            if (pKeep == nullptr) {
                g.mConn = rBlk.mRows;
                return g;
            }
            g.mAllCells = false;
            for (const std::size_t r : *pKeep) {
                g.mConn.insert(g.mConn.end(), rBlk.mRows.begin() + r * rBlk.mWidth,
                               rBlk.mRows.begin() + (r + 1) * rBlk.mWidth);
                g.mCells.push_back(rBlk.mFirst + r);
            }
            return g;
        };
        if (labels.size() > 1) {
            // One block per part and cell block: a rank's cells, as DOLFINx writes them.
            for (const std::int64_t label : labels)
                for (std::size_t k = 0; k < blocks.size(); ++k) {
                    std::vector<std::size_t> keep;
                    for (std::size_t r = 0; r < blocks[k].mCount; ++r)
                        if (parts[k][r] == label)
                            keep.push_back(r);
                    if (!keep.empty())
                        mGroups.push_back(from_rows(blocks[k], &keep));
                }
        } else if (types.size() > 1) {
            for (const Block& blk : blocks) {
                std::vector<std::size_t> all(blk.mCount);
                for (std::size_t r = 0; r < all.size(); ++r)
                    all[r] = r;
                mGroups.push_back(from_rows(blk, &all));
            }
        } else {
            VwGroup g;
            g.mVtkType = blocks.back().mVtk;
            g.mWidth = blocks.back().mWidth;
            for (const Block& blk : blocks)
                g.mConn.insert(g.mConn.end(), blk.mRows.begin(), blk.mRows.end());
            mGroups.push_back(std::move(g));
        }
        mShared = mGroups.size() > 1;
        if (mShared)
            SharePoints();

        std::set<std::string> point_names;
        for (const std::string& name : rMesh.PointDataNames()) {
            if (mShared && (name == kVwGhostType || name == kVwOriginalIds))
                log::warn("VTX (.bp): point data '{}' is replaced by the ghost array of the blocks",
                          name);
            else if (Usable(name, rMesh.PointData(name), mN, "point"))
                point_names.insert(name);
        }
        if (mShared) {
            point_names.insert(kVwGhostType);
            point_names.insert(kVwOriginalIds);
        }
        for (const std::string& name : point_names)
            mArrays.push_back({name, false});
        for (const std::string& name : rMesh.CellDataNames()) {
            std::size_t total = 0;
            bool ok = true;
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
                total += rMesh.CellData(name, b).Size();
            // The concatenation must be (cells[, w]).
            for (std::size_t b = 0; b < rMesh.NumCellBlocks() && ok; ++b) {
                const NDArray& a = rMesh.CellData(name, b);
                const std::size_t rows = rMesh.Cells(b).NumCells();
                ok = Usable(name, a, rows, "cell");
            }
            if (ok && total % mCells == 0)
                mArrays.push_back({name, true});
        }
        mSignature = std::to_string(mN) + "|" + std::to_string(mCells);
        for (const VwGroup& g : mGroups)
            mSignature += "|" + std::to_string(g.mVtkType) + "," + std::to_string(g.mWidth) + "," +
                          std::to_string(g.Rows()) + "," + std::to_string(g.mCells.size());
        for (const VwArray& a : mArrays)
            mSignature += "|" + a.mName + (a.mCell ? "c" : "p");
    }

    const std::string& Signature() const { return mSignature; }

    bool Shared() const { return mShared; }

    // True when `rOther` has this model's points and cells exactly.
    bool SameMesh(const VwModel& rOther) const {
        if (mGeometry != rOther.mGeometry || mGroups.size() != rOther.mGroups.size())
            return false;
        for (std::size_t g = 0; g < mGroups.size(); ++g) {
            const VwGroup& a = mGroups[g];
            const VwGroup& b = rOther.mGroups[g];
            if (a.mVtkType != b.mVtkType || a.mConn != b.mConn || a.mAllPoints != b.mAllPoints ||
                a.mPoints != b.mPoints)
                return false;
        }
        return true;
    }

    std::string Schema() const {
        auto arrays = [&](bool Cell) {
            std::string out;
            for (const VwArray& a : mArrays)
                if (a.mCell == Cell)
                    out += "        <DataArray Name=" + vw_quoteattr(a.mName) + " />\n";
            return out;
        };
        bool any_cell = false;
        for (const VwArray& a : mArrays)
            any_cell = any_cell || a.mCell;
        std::string xml =
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
            "  <UnstructuredGrid>\n"
            "    <Piece NumberOfPoints=\"NumberOfNodes\" NumberOfCells=\"NumberOfEntities\">\n"
            "      <Points>\n"
            "        <DataArray Name=\"geometry\" />\n"
            "      </Points>\n"
            "      <Cells>\n"
            "        <DataArray Name=\"connectivity\" />\n"
            "        <DataArray Name=\"types\" />\n"
            "      </Cells>\n"
            "      <PointData>\n"
            "        <DataArray Name=\"TIME\">step</DataArray>\n" +
            arrays(false) + "      </PointData>\n";
        if (any_cell)
            xml += "      <CellData>\n" + arrays(true) + "      </CellData>\n";
        xml += "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
        return xml;
    }

    void DefineAttributes(adios2::IO& rIo) const {
        rIo.DefineAttribute<std::string>("vtk.xml", Schema());
        if (mShared)
            rIo.DefineAttribute<std::string>(kVtxSharedAttribute, kVwOriginalIds);
    }

    // One step. `WithMesh` false is the step of a `reuse` file after the first: the time and
    // the arrays only.
    void PutStep(adios2::IO& rIo, adios2::Engine& rEngine, const Mesh& rMesh, double Time,
                 bool WithMesh = true) const {
        rEngine.BeginStep();
        Put<double>(rIo, rEngine, "step", &Time, {}, true);
        if (WithMesh) {
            for (int which = 0; which < 2; ++which)
                for (const VwGroup& g : mGroups) {
                    const std::uint32_t value = static_cast<std::uint32_t>(
                        which == 0 ? (g.mAllPoints ? mN : g.mPoints.size()) : g.Rows());
                    const char* name = which == 0 ? "NumberOfNodes" : "NumberOfEntities";
                    if (mShared)
                        Put<std::uint32_t>(rIo, rEngine, name, &value, {1}, false);
                    else
                        PutGlobal(rIo, rEngine, name, value);
                }
            for (const VwGroup& g : mGroups) {
                if (g.mAllPoints) {
                    Put<double>(rIo, rEngine, "geometry", mGeometry.data(), {mN, 3}, false);
                } else {
                    std::vector<double> block;
                    block.reserve(g.mPoints.size() * 3);
                    for (const std::int64_t p : g.mPoints)
                        block.insert(block.end(), mGeometry.begin() + p * 3,
                                     mGeometry.begin() + p * 3 + 3);
                    Put<double>(rIo, rEngine, "geometry", block.data(), {g.mPoints.size(), 3},
                                false);
                }
            }
            for (const VwGroup& g : mGroups)
                Put<std::int64_t>(rIo, rEngine, "connectivity", g.mConn.data(),
                                  {g.Rows(), g.mWidth}, false);
            for (const VwGroup& g : mGroups) {
                const std::uint32_t type = static_cast<std::uint32_t>(g.mVtkType);
                Put<std::uint32_t>(rIo, rEngine, "types", &type, {1}, !mShared);
            }
        }
        for (const VwArray& a : mArrays) {
            if (a.mCell)
                continue;
            const bool ghost_array =
                mShared && (a.mName == kVwGhostType || a.mName == kVwOriginalIds);
            if (!WithMesh && ghost_array)
                continue;
            for (const VwGroup& g : mGroups) {
                if (a.mName == kVwGhostType && mShared) {
                    Put<std::uint8_t>(rIo, rEngine, a.mName, g.mGhost.data(), {g.mPoints.size()},
                                      false);
                } else if (a.mName == kVwOriginalIds && mShared) {
                    Put<std::int64_t>(rIo, rEngine, a.mName, g.mPoints.data(), {g.mPoints.size()},
                                      false);
                } else {
                    const NDArray& arr = rMesh.PointData(a.mName);
                    if (g.mAllPoints) {
                        PutArray(rIo, rEngine, a.mName, arr.Dtype(), arr.Data(), Shape(arr, mN));
                    } else {
                        std::vector<std::size_t> index(g.mPoints.begin(), g.mPoints.end());
                        const std::vector<std::byte> bytes = vw_take_rows(arr, mN, index);
                        PutArray(rIo, rEngine, a.mName, arr.Dtype(), bytes.data(),
                                 Shape(arr, mN, index.size()));
                    }
                }
            }
        }
        for (const VwArray& a : mArrays) {
            if (!a.mCell)
                continue;
            // The blocks' values one after the other, in cell order.
            const NDArray& first = rMesh.CellData(a.mName, 0);
            std::vector<std::byte> bytes;
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                const NDArray& arr = rMesh.CellData(a.mName, b);
                const std::byte* p = static_cast<const std::byte*>(arr.Data());
                bytes.insert(bytes.end(), p, p + arr.Size() * dtype_size(arr.Dtype()));
            }
            const std::size_t row_bytes = bytes.size() / mCells;
            for (const VwGroup& g : mGroups) {
                if (g.mAllCells) {
                    PutArray(rIo, rEngine, a.mName, first.Dtype(), bytes.data(),
                             Shape(first, mCells, mCells));
                } else {
                    std::vector<std::byte> block;
                    block.reserve(g.mCells.size() * row_bytes);
                    for (const std::size_t c : g.mCells)
                        block.insert(block.end(), bytes.begin() + c * row_bytes,
                                     bytes.begin() + (c + 1) * row_bytes);
                    PutArray(rIo, rEngine, a.mName, first.Dtype(), block.data(),
                             Shape(first, mCells, g.mCells.size()));
                }
            }
        }
        rEngine.EndStep();
    }

private:
    // Give every block the points it uses: the first block to use a point owns it, the others
    // hold a ghost copy; the points no cell uses go to the first.
    void SharePoints() {
        std::vector<std::int64_t> owner(mN, -1);
        std::vector<std::vector<std::int64_t>> used(mGroups.size());
        for (std::size_t g = 0; g < mGroups.size(); ++g) {
            std::set<std::int64_t> nodes;
            const VwGroup& grp = mGroups[g];
            for (std::size_t r = 0; r < grp.Rows(); ++r)
                for (std::size_t j = 1; j < grp.mWidth; ++j)
                    nodes.insert(grp.mConn[r * grp.mWidth + j]);
            used[g].assign(nodes.begin(), nodes.end());
            for (const std::int64_t p : used[g])
                if (owner[p] < 0)
                    owner[p] = static_cast<std::int64_t>(g);
        }
        for (std::int64_t& rO : owner)
            if (rO < 0)
                rO = 0;
        for (std::size_t g = 0; g < mGroups.size(); ++g) {
            VwGroup& grp = mGroups[g];
            std::set<std::int64_t> points(used[g].begin(), used[g].end());
            for (std::size_t p = 0; p < mN; ++p)
                if (owner[p] == static_cast<std::int64_t>(g))
                    points.insert(static_cast<std::int64_t>(p));
            grp.mAllPoints = false;
            grp.mPoints.assign(points.begin(), points.end());
            grp.mGhost.clear();
            std::unordered_map<std::int64_t, std::int64_t> local;
            local.reserve(grp.mPoints.size());
            for (std::size_t k = 0; k < grp.mPoints.size(); ++k) {
                local[grp.mPoints[k]] = static_cast<std::int64_t>(k);
                grp.mGhost.push_back(owner[grp.mPoints[k]] != static_cast<std::int64_t>(g) ? 1 : 0);
            }
            for (std::size_t r = 0; r < grp.Rows(); ++r)
                for (std::size_t j = 1; j < grp.mWidth; ++j) {
                    std::int64_t& rV = grp.mConn[r * grp.mWidth + j];
                    rV = local[rV];
                }
        }
    }

    static std::vector<std::size_t> Shape(const NDArray& rArray, std::size_t Rows) {
        return Shape(rArray, Rows, Rows);
    }

    // The block's shape: `Count` rows of the array's `Rows`, with its trailing width.
    static std::vector<std::size_t> Shape(const NDArray& rArray, std::size_t Rows,
                                          std::size_t Count) {
        if (rArray.Shape().size() <= 1)
            return {Count};
        return {Count, rArray.Size() / Rows};
    }

    bool Usable(const std::string& rName, const NDArray& rArray, std::size_t Rows,
                const char* pWhat) const {
        if (vw_reserved(rName)) {
            log::warn("VTX (.bp): {} data '{}' is a reserved name and is not written", pWhat,
                      rName);
            return false;
        }
        if (rArray.Shape().size() > 2 || rArray.Shape().empty() || Rows == 0 ||
            rArray.Size() % Rows != 0) {
            log::warn("VTX (.bp): {} data '{}' is not written", pWhat, rName);
            return false;
        }
        return true;
    }

    // A single value (`Value`) or a rank-local array of `rCount`; a variable written again in
    // the same step is a further block.
    template <class T>
    static void Put(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName,
                    const T* pData, const std::vector<std::size_t>& rCount, bool Value) {
        adios2::Variable<T> var = rIo.InquireVariable<T>(rName);
        if (!var)
            var =
                Value ? rIo.DefineVariable<T>(rName) : rIo.DefineVariable<T>(rName, {}, {}, rCount);
        else if (!Value)
            var.SetSelection({{}, rCount});
        rEngine.Put(var, pData, adios2::Mode::Sync);
    }

    static void PutGlobal(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName,
                          std::uint32_t Value) {
        adios2::Variable<std::uint32_t> var = rIo.InquireVariable<std::uint32_t>(rName);
        if (!var)
            var = rIo.DefineVariable<std::uint32_t>(rName, {1}, {0}, {1});
        rEngine.Put(var, &Value, adios2::Mode::Sync);
    }

    static void PutArray(adios2::IO& rIo, adios2::Engine& rEngine, const std::string& rName,
                         DType Type, const void* pData, const std::vector<std::size_t>& rCount) {
        switch (Type) {
            case DType::Float32:
                return Put<float>(rIo, rEngine, rName, static_cast<const float*>(pData), rCount,
                                  false);
            case DType::Float64:
                return Put<double>(rIo, rEngine, rName, static_cast<const double*>(pData), rCount,
                                   false);
            case DType::Int8:
                return Put<std::int8_t>(rIo, rEngine, rName, static_cast<const std::int8_t*>(pData),
                                        rCount, false);
            case DType::Int16:
                return Put<std::int16_t>(rIo, rEngine, rName,
                                         static_cast<const std::int16_t*>(pData), rCount, false);
            case DType::Int32:
                return Put<std::int32_t>(rIo, rEngine, rName,
                                         static_cast<const std::int32_t*>(pData), rCount, false);
            case DType::Int64:
                return Put<std::int64_t>(rIo, rEngine, rName,
                                         static_cast<const std::int64_t*>(pData), rCount, false);
            case DType::UInt8:
                return Put<std::uint8_t>(rIo, rEngine, rName,
                                         static_cast<const std::uint8_t*>(pData), rCount, false);
            case DType::UInt16:
                return Put<std::uint16_t>(rIo, rEngine, rName,
                                          static_cast<const std::uint16_t*>(pData), rCount, false);
            case DType::UInt32:
                return Put<std::uint32_t>(rIo, rEngine, rName,
                                          static_cast<const std::uint32_t*>(pData), rCount, false);
            case DType::UInt64:
                return Put<std::uint64_t>(rIo, rEngine, rName,
                                          static_cast<const std::uint64_t*>(pData), rCount, false);
        }
    }

    std::size_t mN = 0, mCells = 0;
    bool mShared = false;
    std::vector<double> mGeometry;
    std::vector<VwGroup> mGroups;
    std::vector<VwArray> mArrays;
    std::string mSignature;
};

double vw_time_of(const Mesh& rMesh) {
    if (!rMesh.HasFieldData(kSequenceTimeKey))
        return 0.0;
    const NDArray& a = rMesh.FieldData(kSequenceTimeKey);
    return a.Size() == 0 ? 0.0 : detail::DoubleView(a).Data()[0];
}

}  // namespace

void write_vtx(const std::string& rPath, const Mesh& rMesh) {
    detail::LibraryLock lock;
    const VwModel model(rMesh);
    try {
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("meshioplusplus_vtx_write");
        io.SetEngine("BP5");
        adios2::Engine engine = io.Open(rPath, adios2::Mode::Write);
        model.DefineAttributes(io);
        model.PutStep(io, engine, rMesh, vw_time_of(rMesh));
        engine.Close();
    } catch (const WriteError&) {
        throw;
    } catch (const std::exception& rE) {
        vw_fail(std::string(rE.what()));
    }
}

struct VtxSeriesWriter::Impl {
    std::string mPath;
    VtxMeshPolicy mPolicy = VtxMeshPolicy::Update;
    std::unique_ptr<adios2::ADIOS> mAdios;
    adios2::IO mIo;
    adios2::Engine mEngine;
    std::unique_ptr<VwModel> mModel;
    std::size_t mSteps = 0;
    bool mDone = false;
};

VtxSeriesWriter::VtxSeriesWriter(const std::string& rPath)
    : VtxSeriesWriter(rPath, VtxMeshPolicy::Update) {}

VtxSeriesWriter::VtxSeriesWriter(const std::string& rPath, VtxMeshPolicy Policy)
    : mpImpl(std::make_unique<Impl>()) {
    detail::LibraryLock lock;
    mpImpl->mPath = rPath;
    mpImpl->mPolicy = Policy;
}

VtxSeriesWriter::~VtxSeriesWriter() {
    if (!mpImpl)
        return;
    detail::LibraryLock lock;
    try {
        if (mpImpl->mEngine && !mpImpl->mDone)
            mpImpl->mEngine.Close();
    } catch (...) {
    }
    mpImpl.reset();
}

VtxSeriesWriter::VtxSeriesWriter(VtxSeriesWriter&&) noexcept = default;
VtxSeriesWriter& VtxSeriesWriter::operator=(VtxSeriesWriter&&) noexcept = default;

void VtxSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        vw_fail("series: the writer was moved from");
    detail::LibraryLock lock;
    Impl& s = *mpImpl;
    auto model = std::make_unique<VwModel>(rMesh);
    const bool reuse = s.mPolicy == VtxMeshPolicy::Reuse;
    try {
        if (!s.mModel) {
            s.mAdios = std::make_unique<adios2::ADIOS>();
            s.mIo = s.mAdios->DeclareIO("meshioplusplus_vtx_series");
            s.mIo.SetEngine("BP5");
            s.mEngine = s.mIo.Open(s.mPath, adios2::Mode::Write);
            model->DefineAttributes(s.mIo);
            s.mModel = std::move(model);
        } else if (model->Signature() != s.mModel->Signature()) {
            vw_fail("series: step " + std::to_string(s.mSteps) +
                    "'s points, cells or arrays differ from the first step's; write one file per "
                    "step with '{step}'");
        } else if (reuse && !s.mModel->SameMesh(*model)) {
            vw_fail("series: step " + std::to_string(s.mSteps) +
                    "'s points or cells differ from the first step's, which the reuse policy "
                    "cannot write; use the update policy");
        }
        s.mModel->PutStep(s.mIo, s.mEngine, rMesh, Time, !(reuse && s.mSteps > 0));
        ++s.mSteps;
    } catch (const WriteError&) {
        throw;
    } catch (const std::exception& rE) {
        vw_fail(std::string(rE.what()));
    }
}

std::size_t VtxSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps : 0;
}

void VtxSeriesWriter::Finalize() {
    if (!mpImpl || mpImpl->mDone)
        return;
    detail::LibraryLock lock;
    mpImpl->mDone = true;
    if (mpImpl->mEngine)
        mpImpl->mEngine.Close();
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_ADIOS2
