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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/skin.hpp"

// Project includes (private, not installed)
#include "../detail/closed_facets.hpp"
#include "../detail/crease_edges.hpp"
#include "../detail/step_brep.hpp"
#include "../detail/step_stream.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kStepPrefix = "meshio++: step: ";

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
    write_step(rPath, rMesh, StepWriteOptions{});
}

void write_step(const std::string& rPath, const Mesh& rMesh, const StepWriteOptions& rOptions) {
    detail::write_step_report(rPath, rMesh, rOptions, nullptr);
}

void detail::write_step_report(const std::string& rPath, const Mesh& rMesh,
                               const StepWriteOptions& rOptions, StepBrepResult* pReport,
                               const std::vector<StepFreeformFace>* pFreeform,
                               const std::vector<StepEdgeCurve>* pEdgeCurves, bool BuildBrep) {
    if (rOptions.mFaces != StepFaces::Auto && rOptions.mFaces != StepFaces::Planar &&
        rOptions.mFaces != StepFaces::Fitted && rOptions.mFaces != StepFaces::Facet)
        throw WriteError(std::string(kStepPrefix) +
                         "faces must be 'auto', 'planar', 'fitted' or 'facet'");
    if (std::isnan(rOptions.mTolerance) || std::isinf(rOptions.mTolerance))
        throw WriteError(std::string(kStepPrefix) + "the tolerance must be finite");
    // The facets: the skin of the volume cells when there are any.
    Mesh skin;
    const bool volume = has_skinnable_cells(rMesh);
    if (volume)
        skin = extract_skin(rMesh, /*linearize=*/true);
    const Mesh& surface = volume ? skin : rMesh;
    detail::ClosedFacets shell;
    try {
        shell = detail::closed_facets(surface, kStepPrefix);
    } catch (const std::invalid_argument& rErr) {
        throw WriteError(rErr.what());
    }
    if (shell.mNumSkipped > 0) {
        log::warn("step: {} cell(s) that are not surface facets (lines, vertices) are dropped",
                  shell.mNumSkipped);
        detail::provenance_note("cells-dropped", "cells that are not surface facets");
    }
    if (shell.NumFacets() == 0)
        throw WriteError(std::string(kStepPrefix) + "no surface facets to write");

    // One face per surface: coplanar groups, and the fitted patches.
    const bool fitted_data =
        surface.HasCellData("segment:patch") && surface.HasFieldData("fit:patch_id");
    StepFaces faces = rOptions.mFaces;
    if (faces == StepFaces::Auto)
        faces = fitted_data ? StepFaces::Fitted : StepFaces::Planar;
    if (faces != StepFaces::Facet) {
        detail::StepBrepOptions brep;
        brep.mFitted = faces == StepFaces::Fitted;
        brep.mTolerance = rOptions.mTolerance;
        brep.mFreeform = pFreeform;
        brep.mEdgeCurves = pEdgeCurves;
        brep.mBrep = BuildBrep;
        detail::StepBrepResult r = detail::step_brep_write(rPath, surface, shell, brep);
        if (pReport != nullptr)
            *pReport = std::move(r);
        return;
    }
    const std::vector<double>& xyz = shell.mXyz;
    const std::vector<std::int64_t>& start = shell.mStart;
    const std::vector<std::int64_t>& nodes = shell.mNodes;
    const std::vector<double>& normals = shell.mNormals;
    const std::size_t nf = shell.NumFacets();

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

    detail::StepStream st(rPath, "meshio++ faceted B-rep");
    auto ent = [&](const std::string& rBody) { return st.Ent(rBody); };
    auto ref = [](std::int64_t Id) { return detail::StepStream::Ref(Id); };
    auto triple = [](const double* pV) { return detail::StepStream::Triple(pV); };

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
        const std::int64_t vec = ent("VECTOR(''," + ref(dir) + "," + detail::step_real(len) + ")");
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
    std::vector<std::int64_t> solids;
    for (std::size_t c = 0; c < npieces; ++c) {
        std::string list;
        for (std::size_t f = 0; f < nf; ++f)
            if (piece_of[f] == c)
                list += (list.empty() ? "" : ",") + ref(face_id[f]);
        const std::int64_t shell_id = ent("CLOSED_SHELL('',(" + list + "))");
        solids.push_back(
            ent("MANIFOLD_SOLID_BREP('solid" + std::to_string(c + 1) + "'," + ref(shell_id) + ")"));
    }
    st.Finish(solids);
}

}  // namespace meshioplusplus
