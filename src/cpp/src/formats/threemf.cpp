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
// 3MF packages: the core specification's objects, components, build and base
// materials. See formats/threemf.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/formats/threemf.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/region.hpp"
#include "../detail/open_source.hpp"
#include "../detail/zip_archive.hpp"
#include "pugixml.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kTmfColor = "3mf:color";
constexpr const char* kTmfModelType =
    "http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel";

using TmfMatrix = std::array<double, 12>;  // m00 m01 m02 m10 m11 m12 m20 m21 m22 m30 m31 m32

const TmfMatrix kTmfIdentity = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

/// The 3MF units and their metres.
const std::vector<std::pair<const char*, double>>& tmf_units() {
    static const std::vector<std::pair<const char*, double>> units = {
        {"micron", 1e-6}, {"millimeter", 1e-3}, {"centimeter", 1e-2},
        {"inch", 0.0254}, {"foot", 0.3048},     {"meter", 1.0}};
    return units;
}

std::string_view tmf_local(const char* pName) {
    const std::string_view n(pName);
    const std::size_t colon = n.find(':');
    return colon == std::string_view::npos ? n : n.substr(colon + 1);
}

pugi::xml_node tmf_child(const pugi::xml_node& Node, std::string_view Name) {
    for (const auto c : Node.children())
        if (tmf_local(c.name()) == Name)
            return c;
    return {};
}

double tmf_number(const pugi::xml_node& Node, const char* pAttr) {
    const char* text = Node.attribute(pAttr).value();
    const char* end = nullptr;
    const double v = detail::parse_double(text, end);
    if (end == text || !std::isfinite(v))
        throw ReadError(std::string("meshio++: 3mf: <") + Node.name() + "> has no number '" +
                        pAttr + "'");
    return v;
}

TmfMatrix tmf_transform(const char* pText) {
    if (!pText || !*pText)
        return kTmfIdentity;
    TmfMatrix m{};
    const char* p = pText;
    for (double& v : m) {
        const char* end = nullptr;
        v = detail::parse_double(p, end);
        if (end == p)
            throw ReadError("meshio++: 3mf: a transform needs twelve numbers");
        p = end;
    }
    return m;
}

/// @p A then @p B (a point goes through A first).
TmfMatrix tmf_compose(const TmfMatrix& rA, const TmfMatrix& rB) {
    TmfMatrix out{};
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 3; ++c) {
            double v = r == 3 ? rB[9 + c] : 0.0;
            for (int k = 0; k < 3; ++k)
                v += rA[r * 3 + k] * rB[k * 3 + c];
            out[r * 3 + c] = v;
        }
    return out;
}

std::array<std::uint8_t, 4> tmf_colour(const std::string& rText) {
    std::array<std::uint8_t, 4> rgba = {0, 0, 0, 255};
    if (rText.size() != 7 && rText.size() != 9)
        return rgba;
    for (std::size_t k = 0; k + 1 < rText.size(); k += 2)
        rgba[k / 2] = static_cast<std::uint8_t>(std::stoi(rText.substr(1 + k, 2), nullptr, 16));
    return rgba;
}

struct TmfModel {
    std::map<long long, pugi::xml_node> mObjects;
    std::map<long long, std::vector<std::array<std::uint8_t, 4>>> mMaterials;
};

struct TmfOut {
    std::vector<double> mPoints;
    std::vector<std::int64_t> mTriangles;
    std::vector<std::array<std::uint8_t, 4>> mColours;
    bool mAnyColour = false;
};

void tmf_instance(const TmfModel& rModel, long long Id, const TmfMatrix& rTransform, TmfOut& rOut,
                  int Depth) {
    if (Depth > 64)
        throw ReadError("meshio++: 3mf: components nest deeper than 64 (a cycle?)");
    const auto it = rModel.mObjects.find(Id);
    if (it == rModel.mObjects.end())
        throw ReadError("meshio++: 3mf: the build names object " + std::to_string(Id) +
                        ", which the model does not have");
    const pugi::xml_node object = it->second;
    if (const auto mesh = tmf_child(object, "mesh")) {
        const std::size_t base = rOut.mPoints.size() / 3;
        std::size_t nv = 0;
        for (const auto v : tmf_child(mesh, "vertices").children()) {
            if (tmf_local(v.name()) != "vertex")
                continue;
            const double x = tmf_number(v, "x"), y = tmf_number(v, "y"), z = tmf_number(v, "z");
            const TmfMatrix& m = rTransform;
            rOut.mPoints.push_back(x * m[0] + y * m[3] + z * m[6] + m[9]);
            rOut.mPoints.push_back(x * m[1] + y * m[4] + z * m[7] + m[10]);
            rOut.mPoints.push_back(x * m[2] + y * m[5] + z * m[8] + m[11]);
            ++nv;
        }
        const long long obj_pid = object.attribute("pid").as_llong(-1);
        const long long obj_index = object.attribute("pindex").as_llong(0);
        for (const auto t : tmf_child(mesh, "triangles").children()) {
            if (tmf_local(t.name()) != "triangle")
                continue;
            for (const char* a : {"v1", "v2", "v3"}) {
                const long long v = t.attribute(a).as_llong(-1);
                if (v < 0 || static_cast<std::size_t>(v) >= nv)
                    throw ReadError("meshio++: 3mf: a triangle names vertex " + std::to_string(v) +
                                    " of " + std::to_string(nv));
                rOut.mTriangles.push_back(static_cast<std::int64_t>(base) + v);
            }
            const long long pid = t.attribute("pid").as_llong(obj_pid);
            const long long index = t.attribute("p1").as_llong(t.attribute("pid") ? 0 : obj_index);
            std::array<std::uint8_t, 4> colour = {0, 0, 0, 0};
            const auto mit = rModel.mMaterials.find(pid);
            if (mit != rModel.mMaterials.end() && index >= 0 &&
                static_cast<std::size_t>(index) < mit->second.size()) {
                colour = mit->second[static_cast<std::size_t>(index)];
                rOut.mAnyColour = true;
            }
            rOut.mColours.push_back(colour);
        }
        return;
    }
    if (const auto components = tmf_child(object, "components")) {
        for (const auto c : components.children()) {
            if (tmf_local(c.name()) != "component")
                continue;
            const TmfMatrix local = tmf_transform(c.attribute("transform").value());
            tmf_instance(rModel, c.attribute("objectid").as_llong(-1),
                         tmf_compose(local, rTransform), rOut, Depth + 1);
        }
    }
}

}  // namespace

Mesh read_3mf(const std::string& rPath) {
    const detail::FileSource source =
        detail::open_source(rPath, "meshio++: 3mf: could not open file: " + rPath);
    const detail::ZipReader zip(source.View());
    std::string model_name = "3D/3dmodel.model";
    if (zip.Find("_rels/.rels")) {
        const std::string rels = zip.Read("_rels/.rels");
        pugi::xml_document doc;
        if (doc.load_buffer(rels.data(), rels.size()))
            for (const auto r : doc.document_element().children())
                if (std::string(r.attribute("Type").value()) == kTmfModelType) {
                    model_name = r.attribute("Target").value();
                    if (!model_name.empty() && model_name[0] == '/')
                        model_name.erase(0, 1);
                    break;
                }
    }
    if (!zip.Find(model_name))
        throw ReadError("meshio++: 3mf: the package has no model part '" + model_name + "'");
    const std::string text = zip.Read(model_name);
    pugi::xml_document doc;
    if (!doc.load_buffer(text.data(), text.size()))
        throw ReadError("meshio++: 3mf: the model part is not XML");
    const pugi::xml_node model = doc.document_element();
    if (tmf_local(model.name()) != "model")
        throw ReadError("meshio++: 3mf: the model part has no <model> root");
    const std::string required = model.attribute("requiredextensions").value();
    if (required.find_first_not_of(" \t") != std::string::npos)
        throw ReadError("meshio++: 3mf: the model requires the extensions '" + required +
                        "', which meshio++ does not read (the core specification only)");
    TmfModel m;
    for (const auto r : tmf_child(model, "resources").children()) {
        const std::string_view name = tmf_local(r.name());
        if (name == "object")
            m.mObjects[r.attribute("id").as_llong(-1)] = r;
        else if (name == "basematerials") {
            auto& list = m.mMaterials[r.attribute("id").as_llong(-1)];
            for (const auto b : r.children())
                if (tmf_local(b.name()) == "base")
                    list.push_back(tmf_colour(b.attribute("displaycolor").value()));
        }
    }
    TmfOut out;
    std::vector<std::pair<std::string, std::pair<std::size_t, std::size_t>>> items;
    std::map<std::string, int> seen;
    for (const auto item : tmf_child(model, "build").children()) {
        if (tmf_local(item.name()) != "item")
            continue;
        const long long id = item.attribute("objectid").as_llong(-1);
        const std::size_t first = out.mTriangles.size() / 3;
        tmf_instance(m, id, tmf_transform(item.attribute("transform").value()), out, 0);
        std::string name;
        const auto it = m.mObjects.find(id);
        if (it != m.mObjects.end())
            name = it->second.attribute("name").value();
        if (name.empty())
            name = "object_" + std::to_string(id);
        if (const int n = seen[name]++; n > 0)
            name += "_" + std::to_string(n);
        items.push_back({name, {first, out.mTriangles.size() / 3}});
    }
    Mesh mesh;
    const std::size_t np = out.mPoints.size() / 3, nt = out.mTriangles.size() / 3;
    NDArray pts(DType::Float64, {np, 3});
    std::copy(out.mPoints.begin(), out.mPoints.end(), pts.As<double>());
    mesh.AssignPoints(std::move(pts));
    if (nt > 0) {
        NDArray tri(DType::Int64, {nt, 3});
        std::copy(out.mTriangles.begin(), out.mTriangles.end(), tri.As<std::int64_t>());
        mesh.AddCellBlock("triangle", std::move(tri));
        if (out.mAnyColour) {
            NDArray rgba(DType::UInt8, {nt, 4});
            for (std::size_t t = 0; t < nt; ++t)
                std::memcpy(rgba.As<std::uint8_t>() + t * 4, out.mColours[t].data(), 4);
            mesh.AddCellData(kTmfColor, std::vector<NDArray>{std::move(rgba)});
        }
    }
    for (const auto& [name, span] : items) {
        NDArray ids(DType::Int64, {span.second - span.first});
        for (std::size_t k = span.first; k < span.second; ++k)
            ids.As<std::int64_t>()[k - span.first] = static_cast<std::int64_t>(k);
        mesh.AddRegion(Region(name, RegionKind::Cell, std::move(ids)));
    }
    const std::string unit = model.attribute("unit").value();
    double metres = 1e-3;
    for (const auto& [name, scale] : tmf_units())
        if (unit == name)
            metres = scale;
    detail::set_mesh_length_unit(mesh, metres);
    return mesh;
}

namespace {

std::string tmf_attr_escape(const std::string& rText) {
    std::string out;
    for (const char c : rText) {
        if (c == '&')
            out += "&amp;";
        else if (c == '<')
            out += "&lt;";
        else if (c == '>')
            out += "&gt;";
        else if (c == '"')
            out += "&quot;";
        else
            out += c;
    }
    return out;
}

std::string tmf_g(double V) {
    char buf[40];
    detail::snprintf_c(buf, sizeof buf, "%.17g", V);
    return buf;
}

}  // namespace

void write_3mf(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::array<std::int64_t, 3>> tris;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (std::string(cb.Type()) != "triangle")
            throw WriteError("meshio++: 3mf: a '" + std::string(cb.Type()) +
                             "' block; 3MF meshes are triangles");
        for (std::size_t r = 0; r < cb.NumCells(); ++r)
            tris.push_back({detail::read_int(cb.Conn(), r * 3),
                            detail::read_int(cb.Conn(), r * 3 + 1),
                            detail::read_int(cb.Conn(), r * 3 + 2)});
    }
    // Each triangle to the first cell region holding it, else to "mesh".
    std::vector<int> owner(tris.size(), -1);
    std::vector<std::string> names;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        bool used = false;
        for (std::size_t k = 0; k < r.mEntries.Size(); ++k) {
            const std::int64_t c = detail::read_int(r.mEntries, k);
            if (c >= 0 && static_cast<std::size_t>(c) < tris.size() &&
                owner[static_cast<std::size_t>(c)] < 0) {
                owner[static_cast<std::size_t>(c)] = static_cast<int>(names.size());
                used = true;
            }
        }
        if (used)
            names.push_back(r.mName);
    }
    const int rest = static_cast<int>(names.size());
    if (std::find(owner.begin(), owner.end(), -1) != owner.end() || tris.empty())
        names.push_back("mesh");
    for (int& o : owner)
        if (o < 0)
            o = rest;
    // Base materials from the colour cell data, in order of first appearance.
    std::vector<std::array<std::uint8_t, 4>> colours;
    std::vector<int> colour_of(tris.size(), -1);
    if (rMesh.HasCellData(kTmfColor) &&
        rMesh.CellDataNumBlocks(kTmfColor) == rMesh.NumCellBlocks()) {
        std::size_t t = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& a = rMesh.CellData(kTmfColor, b);
            const std::size_t n = rMesh.Cells(b).NumCells();
            if (a.Size() != n * 4)
                throw WriteError("meshio++: 3mf: 3mf:color needs four uint8 values per triangle");
            for (std::size_t r = 0; r < n; ++r, ++t) {
                std::array<std::uint8_t, 4> c{};
                for (std::size_t k = 0; k < 4; ++k)
                    c[k] = static_cast<std::uint8_t>(detail::read_int(a, r * 4 + k));
                auto it = std::find(colours.begin(), colours.end(), c);
                if (it == colours.end()) {
                    colours.push_back(c);
                    it = colours.end() - 1;
                }
                colour_of[t] = static_cast<int>(it - colours.begin());
            }
        }
    }
    std::string unit = "millimeter";
    const double metres = detail::mesh_length_unit(rMesh);
    for (const auto& [name, scale] : tmf_units())
        if (metres > 0 && std::fabs(metres - scale) <= 1e-12 * scale)
            unit = name;
    std::string x = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<model unit=\"" + unit +
                    "\" xml:lang=\"en-US\" "
                    "xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n"
                    "<metadata name=\"Application\">meshio++</metadata>\n<resources>\n";
    int next_id = 1;
    const int material_id = colours.empty() ? -1 : next_id++;
    if (material_id > 0) {
        x += "<basematerials id=\"" + std::to_string(material_id) + "\">\n";
        char hex[16];
        for (std::size_t k = 0; k < colours.size(); ++k) {
            detail::snprintf_c(hex, sizeof hex, "#%02X%02X%02X%02X", colours[k][0], colours[k][1],
                               colours[k][2], colours[k][3]);
            x += "<base name=\"color_" + std::to_string(k) + "\" displaycolor=\"" + hex + "\"/>\n";
        }
        x += "</basematerials>\n";
    }
    const std::size_t np = rMesh.NumPoints(), pd = rMesh.PointDim();
    std::vector<int> ids;
    for (std::size_t o = 0; o < names.size(); ++o) {
        const int id = next_id++;
        ids.push_back(id);
        x += "<object id=\"" + std::to_string(id) + "\" type=\"model\" name=\"" +
             tmf_attr_escape(names[o]) + "\">\n<mesh>\n<vertices>\n";
        std::vector<std::int64_t> local(np, -1);
        std::vector<std::size_t> order;
        for (std::size_t t = 0; t < tris.size(); ++t)
            if (owner[t] == static_cast<int>(o))
                for (const std::int64_t v : tris[t])
                    if (local[static_cast<std::size_t>(v)] < 0) {
                        local[static_cast<std::size_t>(v)] =
                            static_cast<std::int64_t>(order.size());
                        order.push_back(static_cast<std::size_t>(v));
                    }
        for (const std::size_t v : order) {
            x += "<vertex";
            for (std::size_t c = 0; c < 3; ++c)
                x += std::string(" ") + "xyz"[c] + "=\"" +
                     tmf_g(c < pd ? detail::read_double(rMesh.Points(), v * pd + c) : 0.0) + "\"";
            x += "/>\n";
        }
        x += "</vertices>\n<triangles>\n";
        for (std::size_t t = 0; t < tris.size(); ++t) {
            if (owner[t] != static_cast<int>(o))
                continue;
            x += "<triangle v1=\"" + std::to_string(local[static_cast<std::size_t>(tris[t][0])]) +
                 "\" v2=\"" + std::to_string(local[static_cast<std::size_t>(tris[t][1])]) +
                 "\" v3=\"" + std::to_string(local[static_cast<std::size_t>(tris[t][2])]) + "\"";
            if (material_id > 0)
                x += " pid=\"" + std::to_string(material_id) + "\" p1=\"" +
                     std::to_string(colour_of[t]) + "\"";
            x += "/>\n";
        }
        x += "</triangles>\n</mesh>\n</object>\n";
    }
    x += "</resources>\n<build>\n";
    for (const int id : ids)
        x += "<item objectid=\"" + std::to_string(id) + "\"/>\n";
    x += "</build>\n</model>\n";
    detail::ZipWriter zip;
    zip.Add("[Content_Types].xml",
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Types "
            "xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
            "<Default Extension=\"rels\" "
            "ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
            "<Default Extension=\"model\" "
            "ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/></Types>\n");
    zip.Add("_rels/.rels",
            std::string("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Relationships "
                        "xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                        "<Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"") +
                kTmfModelType + "\"/></Relationships>\n");
    zip.Add("3D/3dmodel.model", x);
    const std::string bytes = zip.Finish();
    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("meshio++: 3mf: could not open file for writing: " + rPath);
    os.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!os)
        throw WriteError("meshio++: 3mf: failed writing " + rPath);
}

}  // namespace meshioplusplus
