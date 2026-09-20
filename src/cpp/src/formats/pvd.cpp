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
// ParaView collection `.pvd`. The index holds no geometry: writing means one
// `.vtu` per step in a sibling directory and an index listing them; reading means
// choosing a step by `timestep`, then a part by `part`, and delegating every entry
// to the reader for its extension (see formats/pindex_common.hpp), which is why
// `.pvd` -> `.pvtu` -> `.vtu` needs no special case.
//
// Anonymous-namespace helpers are prefixed `pvd_`, as the amalgamation
// concatenates every translation unit.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "pindex_common.hpp"
#include "meshioplusplus/formats/pvd.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

// The value of the single-element field-data array @p pKey, if the mesh has one.
bool pvd_field_scalar(const Mesh& rMesh, const char* pKey, double& rOut) {
    if (!rMesh.HasFieldData(pKey))
        return false;
    const NDArray& a = rMesh.FieldData(pKey);
    if (a.Size() != 1)
        return false;
    switch (a.Dtype()) {
        case DType::Float32:
            rOut = static_cast<double>(*a.As<float>());
            return true;
        case DType::Float64:
            rOut = *a.As<double>();
            return true;
        case DType::Int8:
            rOut = static_cast<double>(*a.As<std::int8_t>());
            return true;
        case DType::Int16:
            rOut = static_cast<double>(*a.As<std::int16_t>());
            return true;
        case DType::Int32:
            rOut = static_cast<double>(*a.As<std::int32_t>());
            return true;
        case DType::Int64:
            rOut = static_cast<double>(*a.As<std::int64_t>());
            return true;
        case DType::UInt8:
            rOut = static_cast<double>(*a.As<std::uint8_t>());
            return true;
        case DType::UInt16:
            rOut = static_cast<double>(*a.As<std::uint16_t>());
            return true;
        case DType::UInt32:
            rOut = static_cast<double>(*a.As<std::uint32_t>());
            return true;
        case DType::UInt64:
            rOut = static_cast<double>(*a.As<std::uint64_t>());
            return true;
    }
    return false;
}

struct pvd_entry {
    double mTime = 0.0;
    bool mHasTimestep = false;
    std::int64_t mPart = 0;
    std::size_t mOrder = 0;
    std::string mGroup;
    std::string mName;
    std::string mFile;
};

struct pvd_index {
    std::vector<pvd_entry> mEntries;
    std::vector<double> mTimes;  // distinct, ascending
};

pvd_index pvd_parse(const std::string& rPath) {
    pugi::xml_document doc;
    const pugi::xml_parse_result parsed = doc.load_file(rPath.c_str());
    if (!parsed)
        throw ReadError("meshio++: pvd: could not parse " + rPath + ": " + parsed.description());
    const pugi::xml_node root = doc.child("VTKFile");
    if (!root)
        throw ReadError("meshio++: pvd: expected tag 'VTKFile': " + rPath);
    const std::string type = root.attribute("type").as_string();
    if (type != "Collection")
        throw ReadError("meshio++: pvd: expected type Collection, got '" + type + "': " + rPath);
    const pugi::xml_node coll = root.child("Collection");
    if (!coll)
        throw ReadError("meshio++: pvd: expected tag 'Collection': " + rPath);

    pvd_index index;
    std::size_t order = 0;
    for (const pugi::xml_node ds : coll.children("DataSet")) {
        pvd_entry e;
        try {
            // A missing part is 0. A missing timestep is resolved below: from the
            // entry's own file when it names one, else 0 as ParaView reads it.
            e.mHasTimestep = static_cast<bool>(ds.attribute("timestep"));
            if (e.mHasTimestep)
                e.mTime = detail::stod_c(ds.attribute("timestep").as_string());
            e.mPart = std::stoll(ds.attribute("part").as_string("0"));
        } catch (const std::exception& ex) {
            throw ReadError("meshio++: pvd: bad timestep/part in " + rPath + ": " + ex.what());
        }
        e.mOrder = order++;
        e.mGroup = ds.attribute("group").as_string("");
        e.mName = ds.attribute("name").as_string("");
        e.mFile = ds.attribute("file").as_string("");
        index.mEntries.push_back(std::move(e));
    }

    // The per-file alternative to `timestep=` (VTK's "time in field data"): a
    // `TimeValue` array in the file's own `<FieldData>`, else our `meshio:time`.
    // Resolved here, for every entry, so a summary and a real read cannot group
    // the steps differently. An index whose entries all carry `timestep=` -- what
    // meshio++ and ParaView write -- never opens a piece for this; one that relies
    // on the file's own time opens those files, narrowed to just the two arrays.
    for (pvd_entry& e : index.mEntries) {
        if (e.mHasTimestep || e.mFile.empty())
            continue;
        ReadOptions probe;
        probe.mDataArrays = std::vector<std::string>{"TimeValue", kSequenceTimeKey};
        const Mesh piece = pidx::read_child(pidx::resolve_path(rPath, e.mFile, "pvd", "file"),
                                            probe, /*Wide=*/true);
        if (!pvd_field_scalar(piece, "TimeValue", e.mTime))
            pvd_field_scalar(piece, kSequenceTimeKey, e.mTime);
    }
    for (const pvd_entry& e : index.mEntries)
        index.mTimes.push_back(e.mTime);
    std::sort(index.mTimes.begin(), index.mTimes.end());
    index.mTimes.erase(std::unique(index.mTimes.begin(), index.mTimes.end()), index.mTimes.end());
    return index;
}

std::string pvd_region_name(const pvd_entry& rEntry) {
    if (!rEntry.mName.empty())
        return rEntry.mName;
    const std::string base = "part_" + std::to_string(rEntry.mPart);
    return rEntry.mGroup.empty() ? base : rEntry.mGroup + "/" + base;
}

// The entries of the step at index @p Step, part ascending then document order.
std::vector<const pvd_entry*> pvd_step_entries(const pvd_index& rIndex, std::size_t Step) {
    std::vector<const pvd_entry*> chosen;
    for (const pvd_entry& e : rIndex.mEntries)
        if (e.mTime == rIndex.mTimes[Step])
            chosen.push_back(&e);
    std::stable_sort(chosen.begin(), chosen.end(), [](const pvd_entry* pA, const pvd_entry* pB) {
        return pA->mPart != pB->mPart ? pA->mPart < pB->mPart : pA->mOrder < pB->mOrder;
    });
    return chosen;
}

// The shortest of %.15g / %.16g / %.17g that reads back exactly: locale-independent
// (`snprintf_c` / `stod_c`) and "0.1" rather than "0.10000000000000001".
std::string pvd_format_time(double Time) {
    char buf[40];
    for (int precision = 15; precision <= 17; ++precision) {
        detail::snprintf_c(buf, sizeof(buf), "%.*g", precision, Time);
        if (precision == 17 || detail::stod_c(buf) == Time)
            break;
    }
    return buf;
}

// The step time a plain write records: `meshio:time` when it holds one value.
double pvd_mesh_time(const Mesh& rMesh) {
    double t = 0.0;
    pvd_field_scalar(rMesh, kSequenceTimeKey, t);
    return t;
}

}  // namespace

// --- PvdSeriesWriter -----------------------------------------------------------

struct PvdSeriesWriter::Impl {
    std::string mPath;
    fs::path mDir;
    std::string mStem;
    bool mBinary = true;
    detail::VtkCodec mCodec = detail::VtkCodec::Zlib;
    std::vector<std::pair<double, std::string>> mEntries;  // time, index-relative path
    bool mFinalized = false;

    void FlushIndex() const {
        auto os = detail::make_classic_ofstream(mPath, std::ios::binary);
        if (!os)
            throw WriteError("Could not open file for writing: " + mPath);
        os << "<?xml version=\"1.0\"?>\n";
        os << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        os << detail::provenance_render_xml_comment(detail::SlotTier::Block) << "\n";
        os << "<Collection>\n";
        for (const auto& entry : mEntries)
            os << "<DataSet timestep=\"" << pvd_format_time(entry.first) << "\" part=\"0\" file=\""
               << pidx::escape_attr(entry.second) << "\"/>\n";
        os << "</Collection>\n";
        os << "</VTKFile>\n";
        if (!os)
            throw WriteError("Failed while writing: " + mPath);
    }
};

PvdSeriesWriter::PvdSeriesWriter(const std::string& rPath, bool binary, detail::VtkCodec codec)
    : mImpl(std::make_unique<Impl>()) {
    if (binary && codec != detail::VtkCodec::None)
        detail::vtk_codec_require_write(codec);
    const fs::path index_path(rPath);
    mImpl->mPath = rPath;
    mImpl->mStem = index_path.stem().string();
    mImpl->mDir = index_path.parent_path().empty() ? fs::path(mImpl->mStem)
                                                   : index_path.parent_path() / mImpl->mStem;
    mImpl->mBinary = binary;
    mImpl->mCodec = codec;
    std::error_code ec;
    fs::create_directories(mImpl->mDir, ec);
    if (ec)
        throw WriteError("Could not create directory for .pvd pieces: " + mImpl->mDir.string() +
                         ": " + ec.message());
}

PvdSeriesWriter::~PvdSeriesWriter() = default;
PvdSeriesWriter::PvdSeriesWriter(PvdSeriesWriter&&) noexcept = default;
PvdSeriesWriter& PvdSeriesWriter::operator=(PvdSeriesWriter&&) noexcept = default;

void PvdSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mImpl)
        throw WriteError("meshio++: pvd: Write on a moved-from PvdSeriesWriter");
    if (mImpl->mFinalized)
        throw WriteError("meshio++: pvd: Write after Finalize");
    if (!std::isfinite(Time))
        throw WriteError("meshio++: pvd: a step's time must be finite");
    // Width 4 for every file: the step count is not known up front, and a run
    // whose first files were padded differently from its last is worse than one
    // that is merely wider past 10000 steps.
    std::string idx = std::to_string(mImpl->mEntries.size());
    const std::string name =
        mImpl->mStem + "_" + std::string(idx.size() < 4 ? 4 - idx.size() : 0, '0') + idx + ".vtu";
    write_vtu_codec((mImpl->mDir / name).string(), rMesh, mImpl->mBinary,
                    mImpl->mBinary ? mImpl->mCodec : detail::VtkCodec::None);
    mImpl->mEntries.emplace_back(Time, mImpl->mStem + "/" + name);
    mImpl->FlushIndex();
}

std::size_t PvdSeriesWriter::NumSteps() const noexcept {
    return mImpl ? mImpl->mEntries.size() : 0;
}

void PvdSeriesWriter::Finalize() {
    if (!mImpl || mImpl->mFinalized)
        return;
    if (mImpl->mEntries.empty())
        throw WriteError("meshio++: pvd: a collection needs at least one step");
    mImpl->FlushIndex();
    mImpl->mFinalized = true;
}

// --- write -------------------------------------------------------------------------

void write_pvd(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib) {
    write_pvd_codec(rPath, rMesh, binary, zlib ? detail::VtkCodec::Zlib : detail::VtkCodec::None);
}

void write_pvd_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec) {
    PvdSeriesWriter writer(rPath, binary, codec);
    writer.Write(pvd_mesh_time(rMesh), rMesh);
    writer.Finalize();
}

// --- read --------------------------------------------------------------------------

Mesh read_pvd(const std::string& rPath, const ReadOptions& rOpts) {
    const pvd_index index = pvd_parse(rPath);
    if (index.mEntries.empty())
        return pidx::empty_mesh();

    const std::size_t step = rOpts.ResolveTimeStep(index.mTimes.size());
    const std::vector<const pvd_entry*> chosen = pvd_step_entries(index, step);

    auto one = [&](const pvd_entry& rEntry) {
        Mesh m = pidx::read_child(pidx::resolve_path(rPath, rEntry.mFile, "pvd", "file"), rOpts,
                                  /*Wide=*/true);
        if (rOpts.mGhosts == GhostPolicy::Drop)
            return pidx::drop_ghosts(std::move(m));
        return m;
    };

    Mesh out;
    if (rOpts.mPieceSet) {
        out = one(*chosen[rOpts.ResolvePiece(chosen.size())]);
    } else {
        std::vector<Mesh> pieces;
        std::vector<std::string> names;
        pieces.reserve(chosen.size());
        for (const pvd_entry* e : chosen) {
            pieces.push_back(one(*e));
            names.push_back(pvd_region_name(*e));
        }
        out = pidx::merge_pieces(std::move(pieces), names);
    }
    NDArray time(DType::Float64, {1});
    *time.As<double>() = index.mTimes[step];
    out.AddFieldData(kSequenceTimeKey, std::move(time));
    return out;
}

MeshMetadata read_pvd_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    const pvd_index index = pvd_parse(rPath);
    std::vector<MeshMetadata> parts;
    if (!index.mEntries.empty())
        for (const pvd_entry* e : pvd_step_entries(index, 0))
            parts.push_back(pidx::read_child_metadata(
                pidx::resolve_path(rPath, e->mFile, "pvd", "file"), rOpts, /*Wide=*/true));
    MeshMetadata meta = pidx::aggregate_metadata(parts, "pvd");
    meta.mTimeValues = index.mTimes;
    return meta;
}

}  // namespace meshioplusplus
