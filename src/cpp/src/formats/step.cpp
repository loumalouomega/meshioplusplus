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

// STEP faceted B-rep writer. See formats/step.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/skin.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kStepPrefix = "meshio++: step: ";

// A STEP REAL: `%.17g`, with the decimal point the grammar requires and an
// upper-case exponent ("3" -> "3.", "1e-05" -> "1.E-05").
std::string step_real(double V) {
    char buf[64];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", V);
    std::string s(buf);
    const std::size_t e = s.find('e');
    std::string mant = e == std::string::npos ? s : s.substr(0, e);
    std::string expo = e == std::string::npos ? "" : "E" + s.substr(e + 1);
    if (mant.find('.') == std::string::npos)
        mant += '.';
    return mant + expo;
}

// A STEP string body: printable ASCII, apostrophes doubled.
std::string step_string(const std::string& rText) {
    std::string out;
    for (unsigned char c : rText) {
        if (c == '\'')
            out += "''";
        else if (c >= 0x20 && c < 0x7f)
            out += static_cast<char>(c);
        else
            out += '_';
    }
    return out;
}

std::string step_basename(const std::string& rPath) {
    const std::size_t slash = rPath.find_last_of("/\\");
    return slash == std::string::npos ? rPath : rPath.substr(slash + 1);
}

std::int64_t step_find(const std::vector<std::int64_t>& rSorted, std::int64_t Key) {
    return static_cast<std::int64_t>(std::lower_bound(rSorted.begin(), rSorted.end(), Key) -
                                     rSorted.begin());
}

std::size_t step_uf_find(std::vector<std::size_t>& rParent, std::size_t X) {
    while (rParent[X] != X) {
        rParent[X] = rParent[rParent[X]];
        X = rParent[X];
    }
    return X;
}

void step_unit3(const double* pV, double* pOut) {
    const double len = std::sqrt(pV[0] * pV[0] + pV[1] * pV[1] + pV[2] * pV[2]);
    for (int k = 0; k < 3; ++k)
        pOut[k] = len > 0.0 ? pV[k] / len : 0.0;
}

}  // namespace

void write_step(const std::string& rPath, const Mesh& rMesh) {
    // The facets: the skin of the volume cells when there are any.
    Mesh skin;
    const bool volume = has_skinnable_cells(rMesh);
    if (volume)
        skin = extract_skin(rMesh, /*linearize=*/true);
    const Mesh& surface = volume ? skin : rMesh;
    std::size_t total = 0;
    for (const auto cb : surface.CellRange())
        total += cb.NumCells();
    detail::SurfaceRings in;
    try {
        in = detail::surface_rings(surface, std::vector<std::uint8_t>(total, 1), kStepPrefix);
    } catch (const std::invalid_argument& rErr) {
        throw WriteError(rErr.what());
    }
    if (in.NumFaces() < total) {
        log::warn("step: {} cell(s) that are not surface facets (lines, vertices) are dropped",
                  total - in.NumFaces());
        detail::provenance_note("cells-dropped", "cells that are not surface facets");
    }
    if (in.NumFaces() == 0)
        throw WriteError(std::string(kStepPrefix) + "no surface facets to write");

    const std::vector<double> xyz = detail::surface_xyz(surface);
    for (double v : xyz)
        if (!std::isfinite(v))
            throw WriteError(std::string(kStepPrefix) + "a coordinate is not finite");

    // Fan every non-planar facet into triangles from its first corner.
    std::vector<std::int64_t> start{0}, nodes;
    for (std::size_t f = 0; f < in.NumFaces(); ++f) {
        const std::int64_t* ring = in.mNodes.data() + in.mStart[f];
        const std::size_t size = static_cast<std::size_t>(in.mStart[f + 1] - in.mStart[f]);
        bool planar = true;
        if (size > 3) {
            double n[3];
            detail::ring_unit_normal(xyz.data(), ring, size, n);
            const double* p0 = &xyz[static_cast<std::size_t>(ring[0]) * 3];
            double dev = 0.0, extent = 0.0;
            for (std::size_t k = 1; k < size; ++k) {
                const double* p = &xyz[static_cast<std::size_t>(ring[k]) * 3];
                const double d[3] = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};
                const double h = std::fabs(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
                const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                dev = h > dev ? h : dev;
                extent = l > extent ? l : extent;
            }
            planar = !(dev > 1e-9 * extent);
        }
        if (planar) {
            nodes.insert(nodes.end(), ring, ring + size);
            start.push_back(static_cast<std::int64_t>(nodes.size()));
        } else {
            for (std::size_t k = 1; k + 1 < size; ++k) {
                nodes.insert(nodes.end(), {ring[0], ring[k], ring[k + 1]});
                start.push_back(static_cast<std::int64_t>(nodes.size()));
            }
        }
    }
    const std::size_t nf = start.size() - 1;
    std::vector<double> normals(nf * 3);
    for (std::size_t f = 0; f < nf; ++f) {
        detail::ring_unit_normal(xyz.data(), nodes.data() + start[f],
                                 static_cast<std::size_t>(start[f + 1] - start[f]),
                                 &normals[f * 3]);
        if (normals[f * 3] == 0.0 && normals[f * 3 + 1] == 0.0 && normals[f * 3 + 2] == 0.0)
            throw WriteError(std::string(kStepPrefix) + "facet " + std::to_string(f) +
                             " is degenerate (zero area)");
    }

    // A solid's shell is closed and consistently wound.
    {
        const std::vector<detail::CreaseEdge> bad =
            detail::crease_edges(start, nodes, normals, 180.0);
        std::size_t open = 0, nonmanifold = 0, inconsistent = 0;
        std::string first;
        std::size_t listed = 0;
        for (const detail::CreaseEdge& e : bad) {
            const bool b = e.IsBoundary(), m = e.IsNonManifold(), i = e.mInconsistent;
            if (!b && !m && !i)
                continue;
            open += b ? 1 : 0;
            nonmanifold += m ? 1 : 0;
            inconsistent += i ? 1 : 0;
            if (listed < 5) {
                first += (listed ? ", (" : "(") + std::to_string(e.mLo) + ", " +
                         std::to_string(e.mHi) + ")";
                ++listed;
            }
        }
        if (open + nonmanifold + inconsistent > 0)
            throw WriteError(
                std::string(kStepPrefix) +
                "the surface is not a closed, consistently wound solid: " + std::to_string(open) +
                " open, " + std::to_string(nonmanifold) + " non-manifold and " +
                std::to_string(inconsistent) + " inconsistently wound edge(s), first " + first +
                " (point ids); run repair first");
    }

    // The edges, and the edge-connected pieces (one solid each).
    std::vector<std::pair<std::int64_t, std::int64_t>> edges;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = start[f], s = start[f + 1] - b;
        for (std::int64_t k = 0; k < s; ++k) {
            const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
            const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
            edges.emplace_back(std::min(u, v), std::max(u, v));
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    auto edge_index = [&](std::int64_t U, std::int64_t V) {
        const std::pair<std::int64_t, std::int64_t> key(std::min(U, V), std::max(U, V));
        return static_cast<std::size_t>(std::lower_bound(edges.begin(), edges.end(), key) -
                                        edges.begin());
    };
    std::vector<std::size_t> parent(nf);
    for (std::size_t f = 0; f < nf; ++f)
        parent[f] = f;
    {
        std::vector<std::int64_t> owner(edges.size(), -1);
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t b = start[f], s = start[f + 1] - b;
            for (std::int64_t k = 0; k < s; ++k) {
                const std::size_t e = edge_index(nodes[static_cast<std::size_t>(b + k)],
                                                 nodes[static_cast<std::size_t>(b + (k + 1) % s)]);
                if (owner[e] < 0) {
                    owner[e] = static_cast<std::int64_t>(f);
                    continue;
                }
                const std::size_t ra = step_uf_find(parent, f);
                const std::size_t rb = step_uf_find(parent, static_cast<std::size_t>(owner[e]));
                if (ra != rb)
                    parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
    }
    std::vector<std::size_t> piece_of(nf);
    std::vector<std::size_t> roots;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r = step_uf_find(parent, f);
        if (r == f)
            roots.push_back(f);
        piece_of[f] = static_cast<std::size_t>(std::lower_bound(roots.begin(), roots.end(), r) -
                                               roots.begin());
    }
    const std::size_t npieces = roots.size();

    // A piece wound inwards is written reversed.
    std::vector<double> volume_of(npieces, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = start[f], s = start[f + 1] - b;
        const double* p0 = &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b)]) * 3];
        for (std::int64_t k = 1; k + 1 < s; ++k) {
            const double* p1 =
                &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k)]) * 3];
            const double* p2 =
                &xyz[static_cast<std::size_t>(nodes[static_cast<std::size_t>(b + k + 1)]) * 3];
            volume_of[piece_of[f]] +=
                (p0[0] * (p1[1] * p2[2] - p1[2] * p2[1]) - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0]) +
                 p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])) /
                6.0;
        }
    }
    std::size_t flipped = 0;
    for (double v : volume_of)
        flipped += v < 0.0 ? 1 : 0;
    if (flipped > 0)
        log::warn("step: {} piece(s) of the surface are wound inwards; written reversed", flipped);

    // Points used, ascending.
    std::vector<std::int64_t> used(nodes.begin(), nodes.end());
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());

    std::string out;
    out.reserve(256 * (used.size() + edges.size() + nf) + 4096);
    std::int64_t next = 1;
    auto ent = [&](const std::string& rBody) {
        const std::int64_t id = next++;
        out += '#';
        out += std::to_string(id);
        out += '=';
        out += rBody;
        out += ";\n";
        return id;
    };
    auto ref = [](std::int64_t Id) { return "#" + std::to_string(Id); };
    auto triple = [](const double* pV) {
        return "(" + step_real(pV[0]) + "," + step_real(pV[1]) + "," + step_real(pV[2]) + ")";
    };

    const std::string name = step_string(step_basename(rPath));
    out += "ISO-10303-21;\n/*\n";
    for (std::string line : detail::provenance_lines(detail::SlotTier::Block)) {
        for (std::size_t at = line.find("*/"); at != std::string::npos; at = line.find("*/", at))
            line.replace(at, 2, "* /");
        out += line + "\n";
    }
    out += "*/\nHEADER;\n";
    out += "FILE_DESCRIPTION(('meshio++ faceted B-rep'),'2;1');\n";
    out += "FILE_NAME('" + name + "','1970-01-01T00:00:00',(''),(''),'meshio++','meshio++','');\n";
    out += "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n";

    const std::int64_t app =
        ent("APPLICATION_CONTEXT('core data for automotive mechanical design processes')");
    ent("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000," +
        ref(app) + ")");
    const std::int64_t pctx = ent("PRODUCT_CONTEXT(''," + ref(app) + ",'mechanical')");
    const std::int64_t pdctx =
        ent("PRODUCT_DEFINITION_CONTEXT('part definition'," + ref(app) + ",'design')");
    const std::int64_t prod =
        ent("PRODUCT('" + name + "','" + name + "',''," + "(" + ref(pctx) + "))");
    ent("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(" + ref(prod) + "))");
    const std::int64_t form = ent("PRODUCT_DEFINITION_FORMATION('',''," + ref(prod) + ")");
    const std::int64_t pdef =
        ent("PRODUCT_DEFINITION('design',''," + ref(form) + "," + ref(pdctx) + ")");
    const std::int64_t pshape = ent("PRODUCT_DEFINITION_SHAPE('',''," + ref(pdef) + ")");
    const std::int64_t len_unit = ent("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))");
    const std::int64_t ang_unit = ent("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))");
    const std::int64_t sol_unit = ent("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())");
    const std::int64_t unc =
        ent("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07)," + ref(len_unit) +
            ",'distance_accuracy_value','confusion accuracy')");
    const std::int64_t gctx =
        ent("(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((" + ref(unc) +
            "))GLOBAL_UNIT_ASSIGNED_CONTEXT((" + ref(len_unit) + "," + ref(ang_unit) + "," +
            ref(sol_unit) +
            "))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))");
    const std::int64_t o_pt = ent("CARTESIAN_POINT('',(0.,0.,0.))");
    const std::int64_t o_z = ent("DIRECTION('',(0.,0.,1.))");
    const std::int64_t o_x = ent("DIRECTION('',(1.,0.,0.))");
    const std::int64_t origin =
        ent("AXIS2_PLACEMENT_3D(''," + ref(o_pt) + "," + ref(o_z) + "," + ref(o_x) + ")");

    // Vertices: a point and a vertex per used point id.
    std::vector<std::int64_t> point_id(used.size()), vertex_id(used.size());
    for (std::size_t i = 0; i < used.size(); ++i) {
        point_id[i] =
            ent("CARTESIAN_POINT(''," + triple(&xyz[static_cast<std::size_t>(used[i]) * 3]) + ")");
        vertex_id[i] = ent("VERTEX_POINT(''," + ref(point_id[i]) + ")");
    }
    // Edges: a straight EDGE_CURVE from the lower point id to the higher.
    std::vector<std::int64_t> edge_id(edges.size());
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const double* a = &xyz[static_cast<std::size_t>(edges[e].first) * 3];
        const double* b = &xyz[static_cast<std::size_t>(edges[e].second) * 3];
        const double d[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double len = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        if (!(len > 0.0))
            throw WriteError(std::string(kStepPrefix) + "edge (" + std::to_string(edges[e].first) +
                             ", " + std::to_string(edges[e].second) +
                             ") has zero length (coincident points; run clean first)");
        double u[3];
        step_unit3(d, u);
        const std::size_t ia = static_cast<std::size_t>(step_find(used, edges[e].first));
        const std::size_t ib = static_cast<std::size_t>(step_find(used, edges[e].second));
        const std::int64_t dir = ent("DIRECTION(''," + triple(u) + ")");
        const std::int64_t vec = ent("VECTOR(''," + ref(dir) + "," + step_real(len) + ")");
        const std::int64_t line = ent("LINE(''," + ref(point_id[ia]) + "," + ref(vec) + ")");
        edge_id[e] = ent("EDGE_CURVE(''," + ref(vertex_id[ia]) + "," + ref(vertex_id[ib]) + "," +
                         ref(line) + ",.T.)");
    }
    // Faces: an edge loop on a plane through the first corner.
    std::vector<std::int64_t> face_id(nf);
    std::vector<std::int64_t> ring;
    for (std::size_t f = 0; f < nf; ++f) {
        ring.assign(nodes.begin() + start[f], nodes.begin() + start[f + 1]);
        if (volume_of[piece_of[f]] < 0.0)
            std::reverse(ring.begin(), ring.end());
        const std::size_t s = ring.size();
        std::string loop;
        for (std::size_t k = 0; k < s; ++k) {
            const std::int64_t u = ring[k], v = ring[(k + 1) % s];
            const std::int64_t oe = ent("ORIENTED_EDGE('',*,*," + ref(edge_id[edge_index(u, v)]) +
                                        (u < v ? ",.T.)" : ",.F.)"));
            loop += (k ? "," : "") + ref(oe);
        }
        const std::int64_t el = ent("EDGE_LOOP('',(" + loop + "))");
        const std::int64_t bound = ent("FACE_OUTER_BOUND(''," + ref(el) + ",.T.)");
        double n[3];
        detail::ring_unit_normal(xyz.data(), ring.data(), s, n);
        const double* p0 = &xyz[static_cast<std::size_t>(ring[0]) * 3];
        const double* p1 = &xyz[static_cast<std::size_t>(ring[1]) * 3];
        const double d[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const double h = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
        const double t[3] = {d[0] - h * n[0], d[1] - h * n[1], d[2] - h * n[2]};
        double x[3];
        step_unit3(t, x);
        const std::int64_t nz = ent("DIRECTION(''," + triple(n) + ")");
        const std::int64_t nx = ent("DIRECTION(''," + triple(x) + ")");
        const std::size_t i0 = static_cast<std::size_t>(step_find(used, ring[0]));
        const std::int64_t place =
            ent("AXIS2_PLACEMENT_3D(''," + ref(point_id[i0]) + "," + ref(nz) + "," + ref(nx) + ")");
        const std::int64_t plane = ent("PLANE(''," + ref(place) + ")");
        face_id[f] = ent("ADVANCED_FACE('',(" + ref(bound) + ")," + ref(plane) + ",.T.)");
    }
    // One closed shell and one solid per piece.
    std::string solids;
    for (std::size_t c = 0; c < npieces; ++c) {
        std::string list;
        for (std::size_t f = 0; f < nf; ++f)
            if (piece_of[f] == c)
                list += (list.empty() ? "" : ",") + ref(face_id[f]);
        const std::int64_t shell = ent("CLOSED_SHELL('',(" + list + "))");
        const std::int64_t solid =
            ent("MANIFOLD_SOLID_BREP('solid" + std::to_string(c + 1) + "'," + ref(shell) + ")");
        solids += "," + ref(solid);
    }
    const std::int64_t rep = ent("ADVANCED_BREP_SHAPE_REPRESENTATION('',(" + ref(origin) + solids +
                                 ")," + ref(gctx) + ")");
    ent("SHAPE_DEFINITION_REPRESENTATION(" + ref(pshape) + "," + ref(rep) + ")");
    out += "ENDSEC;\nEND-ISO-10303-21;\n";

    auto os = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);
    os << out;
}

}  // namespace meshioplusplus
