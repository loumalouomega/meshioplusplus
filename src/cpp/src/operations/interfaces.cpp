// Conforming adjacency between selected Cell regions.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kRaPrefix = "meshio++: region_adjacency: ";

struct RaFacet {
    std::vector<std::int64_t> mNodes;
    std::int64_t mCell = -1;
    std::int64_t mFacet = -1;
};

struct RaPairFacet {
    RaFacet mA;
    RaFacet mB;
    std::set<std::int64_t> mCellsA;
    std::set<std::int64_t> mCellsB;
};

struct RaOutputFacet {
    std::vector<std::int64_t> mNodes;
    std::int64_t mRegionA = -1;
    std::int64_t mRegionB = -1;
    std::int64_t mCellA = -1;
    std::int64_t mFacetA = -1;
    std::int64_t mCellB = -1;
    std::int64_t mFacetB = -1;
    std::int64_t mSharedCount = 0;
    double mMeasure = 0.0;
};

struct RaOutputBlock {
    std::string mType;
    std::vector<RaOutputFacet> mFacets;
};

bool ra_matches(const Region& rRegion, const RegionSelector& rSelector) {
    return rRegion.mName == rSelector.mName &&
           (rSelector.mKind < 0 || static_cast<std::int32_t>(rRegion.mKind) == rSelector.mKind) &&
           (rSelector.mDim == kRegionAny || rRegion.mDim == rSelector.mDim) &&
           (rSelector.mTag == kRegionAny || rRegion.mTag == rSelector.mTag);
}

std::size_t ra_find_region(const Mesh& rMesh, const RegionSelector& rSelector) {
    std::size_t found = Mesh::npos;
    std::size_t count = 0;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (ra_matches(rMesh.Region(i), rSelector)) {
            found = i;
            ++count;
        }
    if (count != 1) {
        const std::string state = count == 0 ? "does not exist" : "is ambiguous";
        throw std::invalid_argument(std::string("meshio++: ") + "Region '" + rSelector.mName +
                                    "' " + state + "; pin dim and/or tag");
    }
    return found;
}

void ra_add_facet(std::map<std::vector<std::int64_t>, std::vector<RaFacet>>& rFacets,
                  const std::int64_t* pNodes, std::size_t N, std::int64_t Cell, std::int64_t Facet,
                  std::size_t NumPoints) {
    if (N < 2)
        return;
    std::vector<std::int64_t> nodes(pNodes, pNodes + N);
    for (std::int64_t node : nodes)
        if (node < 0 || static_cast<std::size_t>(node) >= NumPoints)
            return;
    std::vector<std::int64_t> key = nodes;
    std::sort(key.begin(), key.end());
    if (std::adjacent_find(key.begin(), key.end()) != key.end())
        return;
    rFacets[key].push_back(RaFacet{std::move(nodes), Cell, Facet});
}

void ra_collect_facets(const Mesh& rMesh, const std::vector<std::vector<std::size_t>>& rMembership,
                       std::map<std::vector<std::int64_t>, std::vector<RaFacet>>& rFacets) {
    const std::size_t num_points = rMesh.NumPoints();
    std::int64_t base = 0;
    std::vector<std::int64_t> row_nodes;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type_name(cb.Type());
        const std::size_t num_cells = cb.NumCells();
        const bool ragged_polygon = cb.IsRagged() && !cb.IsPolyhedron();
        const CellType type = ragged_polygon ? CellType::Polygon : cell_type_from_name(type_name);
        const int dim = cell_type_dimension(type);

        for (std::size_t row = 0; row < num_cells; ++row) {
            const std::int64_t global_cell = base + static_cast<std::int64_t>(row);
            if (rMembership[static_cast<std::size_t>(global_cell)].empty())
                continue;

            if (cb.IsPolyhedron()) {
                for (std::size_t face = 0; face < cb.NumFaces(row); ++face) {
                    const auto [p_nodes, n_nodes] = cb.Face(row, face);
                    ra_add_facet(rFacets, p_nodes, n_nodes, global_cell,
                                 static_cast<std::int64_t>(face), num_points);
                }
                continue;
            }

            if (ragged_polygon || type_name.rfind("polygon", 0) == 0) {
                const std::size_t n = ragged_polygon ? cb.RowSize(row) : cb.NodesPerCell();
                const std::int64_t* p_nodes = nullptr;
                if (ragged_polygon) {
                    p_nodes = cb.Row(row);
                } else {
                    const NDArray& conn = cb.Conn();
                    row_nodes.resize(n);
                    for (std::size_t i = 0; i < n; ++i)
                        row_nodes[i] = detail::read_int(conn, row * n + i);
                    p_nodes = row_nodes.data();
                }
                for (std::size_t edge = 0; edge < n; ++edge) {
                    const std::int64_t edge_nodes[2] = {p_nodes[edge], p_nodes[(edge + 1) % n]};
                    ra_add_facet(rFacets, edge_nodes, 2, global_cell,
                                 static_cast<std::int64_t>(edge), num_points);
                }
                continue;
            }

            const NDArray& conn = cb.Conn();
            const std::size_t nodes_per_cell = cb.NodesPerCell();
            const std::size_t row_offset = row * nodes_per_cell;
            if (dim == 3) {
                const auto& faces = detail::cell_faces(type);
                for (std::size_t face = 0; face < faces.size(); ++face) {
                    const detail::CellFaceDef& def = faces[face];
                    std::vector<std::int64_t> nodes(def.mNumCorners);
                    for (std::size_t k = 0; k < def.mNumCorners; ++k)
                        nodes[k] = detail::read_int(conn, row_offset + def.mNodes[k]);
                    ra_add_facet(rFacets, nodes.data(), nodes.size(), global_cell,
                                 static_cast<std::int64_t>(face), num_points);
                }
            } else if (dim == 2) {
                const auto& edges = detail::cell_edges(type);
                for (std::size_t edge = 0; edge < edges.size(); ++edge) {
                    const detail::CellEdgeDef& def = edges[edge];
                    const std::int64_t nodes[2] = {
                        detail::read_int(conn, row_offset + def.mNodes[0]),
                        detail::read_int(conn, row_offset + def.mNodes[1]),
                    };
                    ra_add_facet(rFacets, nodes, 2, global_cell, static_cast<std::int64_t>(edge),
                                 num_points);
                }
            }
        }
        base += static_cast<std::int64_t>(num_cells);
    }
}

std::array<double, 3> ra_point(const NDArray& rPoints, std::size_t Dim, std::int64_t Id) {
    std::array<double, 3> point{{0.0, 0.0, 0.0}};
    for (std::size_t d = 0; d < Dim && d < 3; ++d)
        point[d] = detail::read_double(rPoints, static_cast<std::size_t>(Id) * Dim + d);
    return point;
}

bool ra_same_facet_geometry(const Mesh& rMeshA, const RaFacet& rFacetA, const Mesh& rMeshB,
                            const RaFacet& rFacetB) {
    if (rFacetA.mNodes.size() != rFacetB.mNodes.size())
        return false;
    for (std::int64_t node : rFacetA.mNodes) {
        const auto it = std::find(rFacetB.mNodes.begin(), rFacetB.mNodes.end(), node);
        if (it == rFacetB.mNodes.end())
            return false;
        if (ra_point(rMeshA.Points(), rMeshA.PointDim(), node) !=
            ra_point(rMeshB.Points(), rMeshB.PointDim(), node))
            return false;
    }
    return true;
}

double ra_measure(const NDArray& rPoints, std::size_t Dim,
                  const std::vector<std::int64_t>& rNodes) {
    if (rNodes.size() == 2) {
        const auto a = ra_point(rPoints, Dim, rNodes[0]);
        const auto b = ra_point(rPoints, Dim, rNodes[1]);
        const double x = b[0] - a[0];
        const double y = b[1] - a[1];
        const double z = b[2] - a[2];
        return std::sqrt(x * x + y * y + z * z);
    }
    if (rNodes.size() < 3)
        return 0.0;
    const auto origin = ra_point(rPoints, Dim, rNodes[0]);
    double area = 0.0;
    for (std::size_t i = 1; i + 1 < rNodes.size(); ++i) {
        const auto b = ra_point(rPoints, Dim, rNodes[i]);
        const auto c = ra_point(rPoints, Dim, rNodes[i + 1]);
        const std::array<double, 3> u{{b[0] - origin[0], b[1] - origin[1], b[2] - origin[2]}};
        const std::array<double, 3> v{{c[0] - origin[0], c[1] - origin[1], c[2] - origin[2]}};
        const double x = u[1] * v[2] - u[2] * v[1];
        const double y = u[2] * v[0] - u[0] * v[2];
        const double z = u[0] * v[1] - u[1] * v[0];
        area += 0.5 * std::sqrt(x * x + y * y + z * z);
    }
    return area;
}

NDArray ra_int_array(const std::vector<std::int64_t>& rValues) {
    NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), out.As<std::int64_t>());
    return out;
}

NDArray ra_double_array(const std::vector<double>& rValues) {
    NDArray out = NDArray::Uninit(DType::Float64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), out.As<double>());
    return out;
}

std::vector<char> ra_region_mask(const Mesh& rMesh, const RegionSelector& rSelector,
                                 RegionKind Kind, const char* pOperation) {
    std::vector<char> mask(
        static_cast<std::size_t>(detail::total_cells(detail::block_bases(rMesh))), 0);
    const std::string block_prefix = "block:";
    if (rSelector.mName.rfind(block_prefix, 0) == 0) {
        if (rSelector.mKind >= 0 && rSelector.mKind != static_cast<std::int32_t>(Kind))
            throw std::invalid_argument(std::string("meshio++: ") + pOperation +
                                        ": block selectors must select Cell regions");
        std::size_t block = 0;
        try {
            const std::string suffix = rSelector.mName.substr(block_prefix.size());
            std::size_t used = 0;
            block = static_cast<std::size_t>(std::stoull(suffix, &used));
            if (used != suffix.size())
                throw std::invalid_argument("invalid block index");
        } catch (const std::exception&) {
            throw std::invalid_argument(std::string("meshio++: ") + pOperation +
                                        ": invalid block selector '" + rSelector.mName + "'");
        }
        if (block >= rMesh.NumCellBlocks())
            throw std::invalid_argument(std::string("meshio++: ") + pOperation +
                                        ": block selector is out of range");
        const auto bases = detail::block_bases(rMesh);
        for (std::int64_t c = bases[block]; c < bases[block + 1]; ++c)
            mask[static_cast<std::size_t>(c)] = 1;
        return mask;
    }

    RegionSelector selector = rSelector;
    if (selector.mKind < 0)
        selector.mKind = static_cast<std::int32_t>(Kind);
    const std::size_t region_id = ra_find_region(rMesh, selector);
    const Region& region = rMesh.Region(region_id);
    if (region.mKind != Kind)
        throw std::invalid_argument(std::string("meshio++: ") + pOperation + ": selector '" +
                                    selector.mName + "' has the wrong region kind");
    for (std::size_t i = 0; i < region.mEntries.Size(); ++i) {
        const std::int64_t cell = detail::read_int(region.mEntries, i);
        if (cell < 0 || static_cast<std::size_t>(cell) >= mask.size())
            throw std::invalid_argument(std::string("meshio++: ") + pOperation +
                                        ": selected Cell region contains an invalid cell id");
        mask[static_cast<std::size_t>(cell)] = 1;
    }
    return mask;
}

std::vector<RaFacet> ra_boundary_facets(const Mesh& rMesh, const std::vector<char>& rMask) {
    std::vector<std::vector<std::size_t>> membership(rMask.size());
    for (std::size_t c = 0; c < rMask.size(); ++c)
        if (rMask[c])
            membership[c].push_back(0);
    std::map<std::vector<std::int64_t>, std::vector<RaFacet>> indexed;
    ra_collect_facets(rMesh, membership, indexed);
    std::vector<RaFacet> out;
    for (auto& [key, owners] : indexed) {
        (void)key;
        std::sort(owners.begin(), owners.end(), [](const RaFacet& rA, const RaFacet& rB) {
            return std::tie(rA.mCell, rA.mFacet) < std::tie(rB.mCell, rB.mFacet);
        });
        owners.erase(
            std::unique(owners.begin(), owners.end(),
                        [](const RaFacet& rA, const RaFacet& rB) { return rA.mCell == rB.mCell; }),
            owners.end());
        if (owners.size() == 1)
            out.push_back(owners.front());
    }
    std::sort(out.begin(), out.end(), [](const RaFacet& rA, const RaFacet& rB) {
        return std::tie(rA.mCell, rA.mFacet, rA.mNodes) < std::tie(rB.mCell, rB.mFacet, rB.mNodes);
    });
    return out;
}

detail::Vec3 ra_normalize(const detail::Vec3& rV) {
    const double norm = std::sqrt(detail::vec3_norm_sq(rV));
    return norm > 0.0 ? detail::vec3_scale(rV, 1.0 / norm) : detail::Vec3{{0.0, 0.0, 0.0}};
}

detail::Vec3 ra_facet_normal(const Mesh& rMesh, const std::vector<std::int64_t>& rNodes) {
    if (rNodes.size() < 2)
        return {{0.0, 0.0, 0.0}};
    if (rNodes.size() == 2) {
        const auto a = ra_point(rMesh.Points(), rMesh.PointDim(), rNodes[0]);
        const auto b = ra_point(rMesh.Points(), rMesh.PointDim(), rNodes[1]);
        return ra_normalize({{b[1] - a[1], a[0] - b[0], 0.0}});
    }
    const auto a = ra_point(rMesh.Points(), rMesh.PointDim(), rNodes[0]);
    detail::Vec3 normal{{0.0, 0.0, 0.0}};
    for (std::size_t i = 1; i + 1 < rNodes.size(); ++i) {
        const auto b = ra_point(rMesh.Points(), rMesh.PointDim(), rNodes[i]);
        const auto c = ra_point(rMesh.Points(), rMesh.PointDim(), rNodes[i + 1]);
        normal = detail::vec3_add(
            normal, detail::vec3_cross(detail::vec3_sub(b, a), detail::vec3_sub(c, a)));
    }
    return ra_normalize(normal);
}

struct IfProjection {
    detail::Vec3 mPoint{{0.0, 0.0, 0.0}};
    std::array<double, 3> mWeights{{0.0, 0.0, 0.0}};
    std::int64_t mSubfacet = 0;
    double mDistanceSq = std::numeric_limits<double>::infinity();
    detail::Vec3 mNormal{{0.0, 0.0, 0.0}};
    std::int64_t mCell = -1;
    std::int64_t mFacet = -1;
    std::vector<std::int64_t> mNodes;
    bool mFound = false;
};

std::array<double, 3> if_barycentric(const detail::Vec3& rP, const detail::Vec3& rA,
                                     const detail::Vec3& rB, const detail::Vec3& rC) {
    const auto v0 = detail::vec3_sub(rB, rA);
    const auto v1 = detail::vec3_sub(rC, rA);
    const auto v2 = detail::vec3_sub(rP, rA);
    const double d00 = detail::vec3_dot(v0, v0);
    const double d01 = detail::vec3_dot(v0, v1);
    const double d11 = detail::vec3_dot(v1, v1);
    const double d20 = detail::vec3_dot(v2, v0);
    const double d21 = detail::vec3_dot(v2, v1);
    const double den = d00 * d11 - d01 * d01;
    if (!(std::fabs(den) > std::numeric_limits<double>::epsilon()))
        return {{1.0, 0.0, 0.0}};
    const double v = (d11 * d20 - d01 * d21) / den;
    const double w = (d00 * d21 - d01 * d20) / den;
    return {{1.0 - v - w, v, w}};
}

IfProjection if_project_facet(const Mesh& rMesh, const detail::Vec3& rPoint,
                              const RaFacet& rFacet) {
    IfProjection best;
    if (rFacet.mNodes.size() < 2)
        return best;
    if (rFacet.mNodes.size() == 2) {
        const auto a = ra_point(rMesh.Points(), rMesh.PointDim(), rFacet.mNodes[0]);
        const auto b = ra_point(rMesh.Points(), rMesh.PointDim(), rFacet.mNodes[1]);
        const auto ab = detail::vec3_sub(b, a);
        const double len2 = detail::vec3_norm_sq(ab);
        double t = len2 > 0.0 ? detail::vec3_dot(detail::vec3_sub(rPoint, a), ab) / len2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        best.mPoint = detail::vec3_add(a, detail::vec3_scale(ab, t));
        best.mDistanceSq = detail::vec3_norm_sq(detail::vec3_sub(rPoint, best.mPoint));
        best.mWeights = {{1.0 - t, t, 0.0}};
    } else {
        const auto a = ra_point(rMesh.Points(), rMesh.PointDim(), rFacet.mNodes[0]);
        for (std::size_t i = 1; i + 1 < rFacet.mNodes.size(); ++i) {
            const auto b = ra_point(rMesh.Points(), rMesh.PointDim(), rFacet.mNodes[i]);
            const auto c = ra_point(rMesh.Points(), rMesh.PointDim(), rFacet.mNodes[i + 1]);
            const detail::PointTriangleHit hit = detail::closest_point_on_triangle(rPoint, a, b, c);
            if (hit.mDistanceSq < best.mDistanceSq) {
                best.mPoint = hit.mPoint;
                best.mDistanceSq = hit.mDistanceSq;
                const auto ab = detail::vec3_sub(b, a);
                const auto ac = detail::vec3_sub(c, a);
                const auto ap = detail::vec3_sub(hit.mPoint, a);
                switch (hit.mFeature) {
                    case detail::TriangleFeature::VertexA:
                        best.mWeights = {{1.0, 0.0, 0.0}};
                        break;
                    case detail::TriangleFeature::VertexB:
                        best.mWeights = {{0.0, 1.0, 0.0}};
                        break;
                    case detail::TriangleFeature::VertexC:
                        best.mWeights = {{0.0, 0.0, 1.0}};
                        break;
                    case detail::TriangleFeature::EdgeAB: {
                        const double den = detail::vec3_norm_sq(ab);
                        const double t = den > 0.0 ? detail::vec3_dot(ap, ab) / den : 0.0;
                        best.mWeights = {{1.0 - t, t, 0.0}};
                        break;
                    }
                    case detail::TriangleFeature::EdgeBC: {
                        const auto edge = detail::vec3_sub(c, b);
                        const double den = detail::vec3_norm_sq(edge);
                        const double t =
                            den > 0.0
                                ? detail::vec3_dot(detail::vec3_sub(hit.mPoint, b), edge) / den
                                : 0.0;
                        best.mWeights = {{0.0, 1.0 - t, t}};
                        break;
                    }
                    case detail::TriangleFeature::EdgeCA: {
                        const auto edge = detail::vec3_sub(a, c);
                        const double den = detail::vec3_norm_sq(edge);
                        const double t =
                            den > 0.0
                                ? detail::vec3_dot(detail::vec3_sub(hit.mPoint, c), edge) / den
                                : 0.0;
                        best.mWeights = {{t, 0.0, 1.0 - t}};
                        break;
                    }
                    case detail::TriangleFeature::Face:
                        best.mWeights = if_barycentric(hit.mPoint, a, b, c);
                        break;
                }
                best.mSubfacet = static_cast<std::int64_t>(i - 1);
            }
        }
    }
    best.mNormal = ra_facet_normal(rMesh, rFacet.mNodes);
    best.mCell = rFacet.mCell;
    best.mFacet = rFacet.mFacet;
    best.mNodes = rFacet.mNodes;
    best.mFound = std::isfinite(best.mDistanceSq);
    return best;
}

detail::Vec3 if_facet_centroid(const Mesh& rMesh, const RaFacet& rFacet) {
    detail::Vec3 centre{{0.0, 0.0, 0.0}};
    for (std::int64_t node : rFacet.mNodes)
        centre = detail::vec3_add(centre, ra_point(rMesh.Points(), rMesh.PointDim(), node));
    return detail::vec3_scale(centre, 1.0 / static_cast<double>(rFacet.mNodes.size()));
}

double if_mean_edge(const Mesh& rMesh, const std::vector<RaFacet>& rFacets) {
    double sum = 0.0;
    std::size_t count = 0;
    for (const RaFacet& facet : rFacets) {
        for (std::size_t i = 0; i < facet.mNodes.size(); ++i) {
            const auto a = ra_point(rMesh.Points(), rMesh.PointDim(), facet.mNodes[i]);
            const auto b = ra_point(rMesh.Points(), rMesh.PointDim(),
                                    facet.mNodes[(i + 1) % facet.mNodes.size()]);
            sum += std::sqrt(detail::vec3_norm_sq(detail::vec3_sub(b, a)));
            ++count;
        }
    }
    return count == 0 ? 0.0 : sum / static_cast<double>(count);
}

struct IfMatch {
    RaFacet mA;
    RaFacet mB;
    double mGapA = 0.0;
    double mGapB = 0.0;
};

Region if_side_region(const std::string& rName,
                      const std::vector<std::pair<std::int64_t, std::int64_t>>& rEntries) {
    std::vector<std::pair<std::int64_t, std::int64_t>> entries = rEntries;
    std::sort(entries.begin(), entries.end());
    entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
    Region out;
    out.mName = rName;
    out.mKind = RegionKind::Side;
    out.mEntries = NDArray::Uninit(DType::Int64, {entries.size(), std::size_t{2}});
    std::int64_t* data = out.mEntries.As<std::int64_t>();
    for (std::size_t i = 0; i < entries.size(); ++i) {
        data[2 * i] = entries[i].first;
        data[2 * i + 1] = entries[i].second;
    }
    return out;
}

struct IfOutputFacet {
    std::vector<std::int64_t> mNodes;
    std::int64_t mCell = -1;
    std::int64_t mFacet = -1;
    std::int64_t mPartnerCell = -1;
    std::int64_t mPartnerFacet = -1;
    double mGap = 0.0;
    double mMeasure = 0.0;
};

Mesh if_make_facet_mesh(const Mesh& rMaster, const std::vector<IfOutputFacet>& rFacets) {
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMaster.Points()));
    for (const std::string& name : rMaster.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMaster.PointData(name)));
    for (const std::string& name : rMaster.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMaster.FieldData(name)));
    for (std::size_t i = 0; i < rMaster.NumRegions(); ++i)
        if (rMaster.Region(i).mKind == RegionKind::Point)
            out.AddRegion(rMaster.Region(i));

    std::map<std::string, std::vector<IfOutputFacet>> blocks;
    for (const IfOutputFacet& facet : rFacets) {
        const std::string type = facet.mNodes.size() == 2   ? "line"
                                 : facet.mNodes.size() == 3 ? "triangle"
                                 : facet.mNodes.size() == 4 ? "quad"
                                                            : "polygon";
        blocks[type].push_back(facet);
    }
    const std::array<std::string, 4> order{{"line", "triangle", "quad", "polygon"}};
    std::vector<std::vector<std::int64_t>> parent_cells, parent_facets, partner_cells,
        partner_facets;
    std::vector<std::vector<double>> gaps, measures;
    // At most one block per facet type.
    for (auto* pList : {&parent_cells, &parent_facets, &partner_cells, &partner_facets})
        pList->reserve(blocks.size());
    gaps.reserve(blocks.size());
    measures.reserve(blocks.size());
    for (const std::string& type : order) {
        const auto it = blocks.find(type);
        if (it == blocks.end())
            continue;
        const auto& cells = it->second;
        std::vector<std::int64_t> conn, pc, pf, qc, qf;
        std::vector<double> gap, measure;
        const std::size_t width = cells.front().mNodes.size();
        for (auto* pList : {&pc, &pf, &qc, &qf})
            pList->reserve(cells.size());
        gap.reserve(cells.size());
        measure.reserve(cells.size());
        if (type != "polygon")
            conn.reserve(cells.size() * width);
        for (const IfOutputFacet& facet : cells) {
            if (type == "polygon")
                continue;
            conn.insert(conn.end(), facet.mNodes.begin(), facet.mNodes.end());
            pc.push_back(facet.mCell);
            pf.push_back(facet.mFacet);
            qc.push_back(facet.mPartnerCell);
            qf.push_back(facet.mPartnerFacet);
            gap.push_back(facet.mGap);
            measure.push_back(facet.mMeasure);
        }
        if (type == "polygon") {
            std::vector<std::int64_t> flat, offsets{0};
            offsets.reserve(cells.size() + 1);
            for (const IfOutputFacet& facet : cells) {
                flat.insert(flat.end(), facet.mNodes.begin(), facet.mNodes.end());
                offsets.push_back(static_cast<std::int64_t>(flat.size()));
                pc.push_back(facet.mCell);
                pf.push_back(facet.mFacet);
                qc.push_back(facet.mPartnerCell);
                qf.push_back(facet.mPartnerFacet);
                gap.push_back(facet.mGap);
                measure.push_back(facet.mMeasure);
            }
            out.AddPolygonBlock("polygon", std::move(flat), std::move(offsets));
        } else {
            NDArray connectivity = NDArray::Uninit(DType::Int64, {cells.size(), width});
            std::copy(conn.begin(), conn.end(), connectivity.As<std::int64_t>());
            out.AddCellBlock(type, std::move(connectivity));
        }
        parent_cells.push_back(std::move(pc));
        parent_facets.push_back(std::move(pf));
        partner_cells.push_back(std::move(qc));
        partner_facets.push_back(std::move(qf));
        gaps.push_back(std::move(gap));
        measures.push_back(std::move(measure));
    }
    if (!parent_cells.empty()) {
        auto add_i64 = [&](const char* pName,
                           const std::vector<std::vector<std::int64_t>>& rValues) {
            std::vector<NDArray> arrays;
            arrays.reserve(rValues.size());
            for (const auto& values : rValues)
                arrays.push_back(ra_int_array(values));
            out.AddCellData(pName, std::move(arrays));
        };
        auto add_f64 = [&](const char* pName, const std::vector<std::vector<double>>& rValues) {
            std::vector<NDArray> arrays;
            arrays.reserve(rValues.size());
            for (const auto& values : rValues)
                arrays.push_back(ra_double_array(values));
            out.AddCellData(pName, std::move(arrays));
        };
        add_i64("interface:parent_cell", parent_cells);
        add_i64("interface:parent_facet", parent_facets);
        add_i64("interface:partner_cell", partner_cells);
        add_i64("interface:partner_facet", partner_facets);
        add_f64("interface:gap", gaps);
        add_f64("interface:measure", measures);
    }
    return out;
}

}  // namespace

Mesh region_adjacency(const Mesh& rMesh, const std::vector<RegionSelector>& rRegions) {
    std::vector<std::size_t> region_indices;
    std::vector<std::string> region_names;
    bool block_groups = false;
    if (rRegions.empty()) {
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
            if (rMesh.Region(i).mKind == RegionKind::Cell)
                region_indices.push_back(i);
        if (region_indices.size() < 2 && rMesh.NumCellBlocks() >= 2) {
            region_indices.clear();
            block_groups = true;
            for (std::size_t i = 0; i < rMesh.NumCellBlocks(); ++i)
                region_names.push_back("block:" + std::to_string(i));
        }
    } else {
        for (const RegionSelector& selector : rRegions) {
            RegionSelector cell_selector = selector;
            if (cell_selector.mKind == -1)
                cell_selector.mKind = static_cast<std::int32_t>(RegionKind::Cell);
            const std::size_t index = ra_find_region(rMesh, cell_selector);
            if (rMesh.Region(index).mKind != RegionKind::Cell)
                throw std::invalid_argument(std::string(kRaPrefix) +
                                            "selectors must name Cell regions");
            if (std::find(region_indices.begin(), region_indices.end(), index) !=
                region_indices.end())
                throw std::invalid_argument(std::string(kRaPrefix) +
                                            "the same region was selected twice");
            region_indices.push_back(index);
        }
    }
    if (region_indices.size() < 2)
        if (!block_groups || region_names.size() < 2)
            throw std::invalid_argument(std::string(kRaPrefix) +
                                        "select at least two Cell regions or use at least two "
                                        "cell blocks");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t num_cells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<std::vector<std::size_t>> membership(num_cells);
    if (block_groups) {
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            for (std::int64_t cell = bases[b]; cell < bases[b + 1]; ++cell)
                membership[static_cast<std::size_t>(cell)].push_back(b);
    } else {
        region_names.reserve(region_indices.size());
        for (std::size_t i = 0; i < region_indices.size(); ++i) {
            const Region& region = rMesh.Region(region_indices[i]);
            region_names.push_back(region.mName);
            for (std::size_t e = 0; e < region.mEntries.Size(); ++e) {
                const std::int64_t cell = detail::read_int(region.mEntries, e);
                if (cell >= 0 && static_cast<std::size_t>(cell) < num_cells)
                    membership[static_cast<std::size_t>(cell)].push_back(i);
            }
        }
    }

    std::map<std::vector<std::int64_t>, std::vector<RaFacet>> facets;
    ra_collect_facets(rMesh, membership, facets);

    using PairKey = std::tuple<std::size_t, std::size_t, std::vector<std::int64_t>>;
    std::map<PairKey, RaPairFacet> pair_facets;
    for (const auto& [key, owners] : facets) {
        for (std::size_t i = 0; i < owners.size(); ++i) {
            for (std::size_t j = i + 1; j < owners.size(); ++j) {
                if (owners[i].mCell == owners[j].mCell)
                    continue;
                const auto& mem_i = membership[static_cast<std::size_t>(owners[i].mCell)];
                const auto& mem_j = membership[static_cast<std::size_t>(owners[j].mCell)];
                for (std::size_t a : mem_i)
                    for (std::size_t b : mem_j) {
                        if (a == b)
                            continue;
                        const bool forward = a < b;
                        const std::size_t region_a = forward ? a : b;
                        const std::size_t region_b = forward ? b : a;
                        RaPairFacet& pair = pair_facets[{region_a, region_b, key}];
                        const RaFacet& facet_a = forward ? owners[i] : owners[j];
                        const RaFacet& facet_b = forward ? owners[j] : owners[i];
                        if (pair.mCellsA.empty() || std::tie(facet_a.mCell, facet_a.mFacet) <
                                                        std::tie(pair.mA.mCell, pair.mA.mFacet))
                            pair.mA = facet_a;
                        if (pair.mCellsB.empty() || std::tie(facet_b.mCell, facet_b.mFacet) <
                                                        std::tie(pair.mB.mCell, pair.mB.mFacet))
                            pair.mB = facet_b;
                        pair.mCellsA.insert(facet_a.mCell);
                        pair.mCellsB.insert(facet_b.mCell);
                    }
            }
        }
    }

    std::map<std::string, RaOutputBlock> output_blocks;
    for (const auto& [pair_key, pair] : pair_facets) {
        const auto [region_a, region_b, key] = pair_key;
        (void)key;
        const std::size_t node_count = pair.mA.mNodes.size();
        const std::string type = node_count == 2   ? "line"
                                 : node_count == 3 ? "triangle"
                                 : node_count == 4 ? "quad"
                                                   : "polygon";
        RaOutputBlock& block = output_blocks[type];
        block.mType = type;
        block.mFacets.push_back(RaOutputFacet{
            pair.mA.mNodes, static_cast<std::int64_t>(region_a),
            static_cast<std::int64_t>(region_b), pair.mA.mCell, pair.mA.mFacet, pair.mB.mCell,
            pair.mB.mFacet, static_cast<std::int64_t>(pair.mCellsA.size() + pair.mCellsB.size()),
            ra_measure(rMesh.Points(), rMesh.PointDim(), pair.mA.mNodes)});
    }

    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mKind == RegionKind::Point)
            out.AddRegion(rMesh.Region(i));

    const std::array<std::string, 4> order{{"line", "triangle", "quad", "polygon"}};
    std::vector<std::vector<std::int64_t>> region_a_blocks, region_b_blocks, parent_a_blocks,
        facet_a_blocks, parent_b_blocks, facet_b_blocks, shared_blocks;
    std::vector<std::vector<double>> measure_blocks;
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> output_regions;
    std::int64_t output_base = 0;
    std::vector<std::int64_t> polygon_flat, polygon_rows{0};
    for (const std::string& type : order) {
        auto it = output_blocks.find(type);
        if (it == output_blocks.end())
            continue;
        RaOutputBlock& block = it->second;
        std::vector<std::int64_t> conn;
        const std::size_t width = block.mFacets.empty() ? 0 : block.mFacets.front().mNodes.size();
        std::vector<std::int64_t> ids_a, ids_b, cells_a, facets_a, cells_b, facets_b, shared;
        std::vector<double> measure;
        for (const RaOutputFacet& facet : block.mFacets) {
            if (type == "polygon") {
                polygon_flat.insert(polygon_flat.end(), facet.mNodes.begin(), facet.mNodes.end());
                polygon_rows.push_back(static_cast<std::int64_t>(polygon_flat.size()));
            } else {
                conn.insert(conn.end(), facet.mNodes.begin(), facet.mNodes.end());
            }
            ids_a.push_back(facet.mRegionA);
            ids_b.push_back(facet.mRegionB);
            cells_a.push_back(facet.mCellA);
            facets_a.push_back(facet.mFacetA);
            cells_b.push_back(facet.mCellB);
            facets_b.push_back(facet.mFacetB);
            shared.push_back(facet.mSharedCount);
            measure.push_back(facet.mMeasure);
        }
        if (type == "polygon") {
            // Polygons may have different arities; row offsets and the flat
            // connectivity move directly into the backend's CSR block.
            std::vector<std::int64_t> local_flat;
            std::vector<std::int64_t> local_rows{0};
            for (const RaOutputFacet& facet : block.mFacets) {
                local_flat.insert(local_flat.end(), facet.mNodes.begin(), facet.mNodes.end());
                local_rows.push_back(static_cast<std::int64_t>(local_flat.size()));
            }
            out.AddPolygonBlock("polygon", std::move(local_flat), std::move(local_rows));
        } else {
            NDArray array = NDArray::Uninit(DType::Int64, {block.mFacets.size(), width});
            std::copy(conn.begin(), conn.end(), array.As<std::int64_t>());
            out.AddCellBlock(type, std::move(array));
        }
        for (const RaOutputFacet& facet : block.mFacets) {
            const std::string pair_name = "adjacency:" + std::to_string(facet.mRegionA) + ":" +
                                          std::to_string(facet.mRegionB);
            auto region_it = std::find_if(output_regions.begin(), output_regions.end(),
                                          [&](const auto& r) { return r.first == pair_name; });
            if (region_it == output_regions.end())
                output_regions.push_back({pair_name, {output_base}});
            else
                region_it->second.push_back(output_base);
            ++output_base;
        }
        region_a_blocks.push_back(std::move(ids_a));
        region_b_blocks.push_back(std::move(ids_b));
        parent_a_blocks.push_back(std::move(cells_a));
        facet_a_blocks.push_back(std::move(facets_a));
        parent_b_blocks.push_back(std::move(cells_b));
        facet_b_blocks.push_back(std::move(facets_b));
        shared_blocks.push_back(std::move(shared));
        measure_blocks.push_back(std::move(measure));
    }

    auto add_int_data = [&](const std::string& rName,
                            const std::vector<std::vector<std::int64_t>>& rBlocks) {
        std::vector<NDArray> arrays;
        arrays.reserve(rBlocks.size());
        for (const auto& block : rBlocks)
            arrays.push_back(ra_int_array(block));
        out.AddCellData(rName, std::move(arrays));
    };
    add_int_data("interface:region_a", region_a_blocks);
    add_int_data("interface:region_b", region_b_blocks);
    add_int_data("interface:parent_cell_a", parent_a_blocks);
    add_int_data("interface:parent_facet_a", facet_a_blocks);
    add_int_data("interface:parent_cell_b", parent_b_blocks);
    add_int_data("interface:parent_facet_b", facet_b_blocks);
    add_int_data("interface:shared_count", shared_blocks);
    std::vector<NDArray> measures;
    measures.reserve(measure_blocks.size());
    for (const auto& block : measure_blocks)
        measures.push_back(ra_double_array(block));
    out.AddCellData("interface:measure", std::move(measures));

    for (auto& [name, cells] : output_regions) {
        Region region;
        region.mName = std::move(name);
        region.mKind = RegionKind::Cell;
        region.mDim = -1;
        region.mTag = -1;
        region.mEntries = ra_int_array(cells);
        out.AddRegion(std::move(region));
    }
    if (rMesh.NumRegions() > 0) {
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
            if (rMesh.Region(i).mKind != RegionKind::Point) {
                log::warn(
                    "region_adjacency: source Cell/Side regions are represented by the "
                    "output adjacency regions; original entries are not copied");
                break;
            }
    }
    return out;
}

FindInterfaceResult find_interface(const Mesh& rMesh, const RegionSelector& rRegionA,
                                   const RegionSelector& rRegionB,
                                   const FindInterfaceOptions& rOptions) {
    return find_interface(rMesh, rRegionA, rMesh, rRegionB, rOptions);
}

FindInterfaceResult find_interface(const Mesh& rMeshA, const RegionSelector& rRegionA,
                                   const Mesh& rMeshB, const RegionSelector& rRegionB,
                                   const FindInterfaceOptions& rOptions) {
    if (rOptions.mMode != InterfaceMode::Conforming && rOptions.mMode != InterfaceMode::Proximity)
        throw std::invalid_argument(
            "meshio++: find_interface: mode must be conforming or proximity");
    if (rOptions.mMaster != InterfaceMaster::A && rOptions.mMaster != InterfaceMaster::B)
        throw std::invalid_argument("meshio++: find_interface: master must be A or B");
    if (!(rOptions.mGapTolerance >= 0.0) || !std::isfinite(rOptions.mGapTolerance))
        throw std::invalid_argument(
            "meshio++: find_interface: gap_tolerance must be finite and non-negative");
    if (!(rOptions.mOverlapTolerance >= 0.0) || !std::isfinite(rOptions.mOverlapTolerance))
        throw std::invalid_argument(
            "meshio++: find_interface: overlap_tolerance must be finite and non-negative");
    if (!(rOptions.mAngleTolerance >= 0.0 && rOptions.mAngleTolerance <= 180.0) ||
        !std::isfinite(rOptions.mAngleTolerance))
        throw std::invalid_argument(
            "meshio++: find_interface: angle_tolerance must be in [0, 180]");

    const auto mask_a = ra_region_mask(rMeshA, rRegionA, RegionKind::Cell, "find_interface");
    const auto mask_b = ra_region_mask(rMeshB, rRegionB, RegionKind::Cell, "find_interface");
    const std::vector<RaFacet> facets_a = ra_boundary_facets(rMeshA, mask_a);
    const std::vector<RaFacet> facets_b = ra_boundary_facets(rMeshB, mask_b);
    if (facets_a.empty() || facets_b.empty())
        throw std::invalid_argument(
            "meshio++: find_interface: both selected parts must have boundary facets");

    std::vector<IfMatch> matches;
    std::set<std::pair<std::int64_t, std::int64_t>> matched_a, matched_b;
    if (rOptions.mMode == InterfaceMode::Conforming) {
        std::map<std::vector<std::int64_t>, std::vector<RaFacet>> map_a, map_b;
        for (const RaFacet& facet : facets_a) {
            auto key = facet.mNodes;
            std::sort(key.begin(), key.end());
            map_a[std::move(key)].push_back(facet);
        }
        for (const RaFacet& facet : facets_b) {
            auto key = facet.mNodes;
            std::sort(key.begin(), key.end());
            map_b[std::move(key)].push_back(facet);
        }
        for (const auto& [key, owners_a] : map_a) {
            const auto it = map_b.find(key);
            if (it == map_b.end())
                continue;
            for (const RaFacet& a : owners_a)
                for (const RaFacet& b : it->second) {
                    if (&rMeshA == &rMeshB && a.mCell == b.mCell)
                        continue;
                    if (!ra_same_facet_geometry(rMeshA, a, rMeshB, b))
                        continue;
                    matches.push_back({a, b, 0.0, 0.0});
                    matched_a.insert({a.mCell, a.mFacet});
                    matched_b.insert({b.mCell, b.mFacet});
                }
        }
    } else {
        double tolerance = rOptions.mGapTolerance;
        if (tolerance == 0.0) {
            const double edge_a = if_mean_edge(rMeshA, facets_a);
            const double edge_b = if_mean_edge(rMeshB, facets_b);
            tolerance = 0.01 * (edge_a + edge_b) * 0.5;
        }
        tolerance += rOptions.mOverlapTolerance;
        const double max_distance_sq = tolerance * tolerance;
        const double normal_cos = std::cos(rOptions.mAngleTolerance * (std::acos(-1.0) / 180.0));
        // The search runs from the master side, so every master facet that has a partner within
        // tolerance is reported: an A-driven search would only return the B facets A happens to
        // pick, which is an incomplete (and duplicated) interface when B is the finer mesh.
        const bool swapped = rOptions.mMaster == InterfaceMaster::B;
        const Mesh& rSearchMesh = swapped ? rMeshB : rMeshA;
        const Mesh& rTargetMesh = swapped ? rMeshA : rMeshB;
        const std::vector<RaFacet>& search_facets = swapped ? facets_b : facets_a;
        const std::vector<RaFacet>& target_facets = swapped ? facets_a : facets_b;
        for (const RaFacet& s : search_facets) {
            const detail::Vec3 centre = if_facet_centroid(rSearchMesh, s);
            const detail::Vec3 normal_s = ra_facet_normal(rSearchMesh, s.mNodes);
            if (!(detail::vec3_norm_sq(normal_s) > 0.0))
                continue;
            IfProjection best;
            for (const RaFacet& t : target_facets) {
                const IfProjection projection = if_project_facet(rTargetMesh, centre, t);
                if (!projection.mFound || projection.mDistanceSq > max_distance_sq ||
                    !(detail::vec3_norm_sq(projection.mNormal) > 0.0))
                    continue;
                if (detail::vec3_dot(normal_s, projection.mNormal) > -normal_cos)
                    continue;
                if (!best.mFound || projection.mDistanceSq < best.mDistanceSq ||
                    (projection.mDistanceSq == best.mDistanceSq &&
                     std::tie(projection.mCell, projection.mFacet) <
                         std::tie(best.mCell, best.mFacet)))
                    best = projection;
            }
            if (!best.mFound)
                continue;
            bool all_corners_near = true;
            for (std::int64_t node : s.mNodes) {
                const detail::Vec3 point =
                    ra_point(rSearchMesh.Points(), rSearchMesh.PointDim(), node);
                double nearest = std::numeric_limits<double>::infinity();
                for (const RaFacet& t : target_facets)
                    nearest =
                        std::min(nearest, if_project_facet(rTargetMesh, point, t).mDistanceSq);
                if (nearest > max_distance_sq) {
                    all_corners_near = false;
                    break;
                }
            }
            if (!all_corners_near)
                continue;
            const detail::Vec3 delta = detail::vec3_sub(centre, best.mPoint);
            const double gap_s = detail::vec3_dot(delta, normal_s);
            const double gap_t = detail::vec3_dot(delta, best.mNormal);
            RaFacet t{best.mNodes, best.mCell, best.mFacet};
            if (swapped)
                matches.push_back({std::move(t), s, gap_t, gap_s});
            else
                matches.push_back({s, std::move(t), gap_s, gap_t});
            matched_a.insert(swapped ? std::make_pair(best.mCell, best.mFacet)
                                     : std::make_pair(s.mCell, s.mFacet));
            matched_b.insert(swapped ? std::make_pair(s.mCell, s.mFacet)
                                     : std::make_pair(best.mCell, best.mFacet));
        }
    }

    std::sort(matches.begin(), matches.end(), [](const IfMatch& rA, const IfMatch& rB) {
        return std::tie(rA.mA.mCell, rA.mA.mFacet, rA.mB.mCell, rA.mB.mFacet) <
               std::tie(rB.mA.mCell, rB.mA.mFacet, rB.mB.mCell, rB.mB.mFacet);
    });
    std::vector<IfOutputFacet> output_facets;
    output_facets.reserve(matches.size());
    std::vector<std::pair<std::int64_t, std::int64_t>> side_a, side_b;
    InterfaceReport report;
    report.mNumPairs = static_cast<std::int64_t>(matches.size());
    report.mUnmatchedA = static_cast<std::int64_t>(facets_a.size() - matched_a.size());
    report.mUnmatchedB = static_cast<std::int64_t>(facets_b.size() - matched_b.size());
    for (const IfMatch& match : matches) {
        side_a.push_back({match.mA.mCell, match.mA.mFacet});
        side_b.push_back({match.mB.mCell, match.mB.mFacet});
        const RaFacet& master = rOptions.mMaster == InterfaceMaster::A ? match.mA : match.mB;
        const RaFacet& partner = rOptions.mMaster == InterfaceMaster::A ? match.mB : match.mA;
        const Mesh& master_mesh = rOptions.mMaster == InterfaceMaster::A ? rMeshA : rMeshB;
        const double measure =
            ra_measure(master_mesh.Points(), master_mesh.PointDim(), master.mNodes);
        const double gap = rOptions.mMode == InterfaceMode::Conforming ? 0.0
                           : rOptions.mMaster == InterfaceMaster::A    ? match.mGapA
                                                                       : match.mGapB;
        output_facets.push_back({master.mNodes, master.mCell, master.mFacet, partner.mCell,
                                 partner.mFacet, gap, measure});
        report.mArea += measure;
        report.mMaxGap = std::max(report.mMaxGap, std::fabs(gap));
    }
    const Mesh& master_mesh = rOptions.mMaster == InterfaceMaster::A ? rMeshA : rMeshB;
    FindInterfaceResult result;
    result.mMesh = if_make_facet_mesh(master_mesh, output_facets);
    result.mSideA = if_side_region("interface:side_a", side_a);
    result.mSideB = if_side_region("interface:side_b", side_b);
    result.mReport = report;
    return result;
}

namespace {

std::vector<std::int64_t> if_self_cell_corners(const Mesh::CellView& rBlock, std::size_t Row) {
    std::vector<std::int64_t> nodes;
    const std::string type_name(rBlock.Type());
    if (rBlock.IsRagged()) {
        if (rBlock.IsPolyhedron())
            return nodes;
        const std::int64_t* row = rBlock.Row(Row);
        nodes.assign(row, row + rBlock.RowSize(Row));
        return nodes;
    }
    const CellType type = cell_type_from_name(type_name);
    const std::size_t width = rBlock.NodesPerCell();
    const int corner_count = detail::cell_corner_count(type);
    if (corner_count <= 0)
        return nodes;
    const std::size_t corners = static_cast<std::size_t>(corner_count);
    const NDArray& conn = rBlock.Conn();
    nodes.reserve(corners);
    for (std::size_t i = 0; i < corners; ++i)
        nodes.push_back(detail::read_int(conn, Row * width + i));
    (void)type;
    return nodes;
}

std::vector<RaFacet> if_contact_facets(const Mesh& rMesh, const std::vector<char>& rMask) {
    const auto bases = detail::block_bases(rMesh);
    std::vector<char> boundary_mask = rMask;
    std::vector<RaFacet> facets;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type_name(cb.Type());
        const CellType type =
            cb.IsPolyhedron() ? CellType::Polyhedron : cell_type_from_name(type_name);
        if (cell_type_dimension(type) != 2 || rMesh.PointDim() < 3)
            continue;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t global = bases[b] + static_cast<std::int64_t>(c);
            if (!rMask[static_cast<std::size_t>(global)])
                continue;
            std::vector<std::int64_t> corners = if_self_cell_corners(cb, c);
            if (!corners.empty())
                facets.push_back({std::move(corners), global, 0});
            boundary_mask[static_cast<std::size_t>(global)] = 0;
        }
    }
    std::vector<RaFacet> boundary = ra_boundary_facets(rMesh, boundary_mask);
    facets.insert(facets.end(), boundary.begin(), boundary.end());
    std::sort(facets.begin(), facets.end(), [](const RaFacet& rA, const RaFacet& rB) {
        return std::tie(rA.mCell, rA.mFacet, rA.mNodes) < std::tie(rB.mCell, rB.mFacet, rB.mNodes);
    });
    return facets;
}

NDArray if_f64_matrix(const std::vector<double>& rValues, std::size_t Rows, std::size_t Cols) {
    NDArray out = NDArray::Uninit(DType::Float64, {Rows, Cols});
    if (!rValues.empty())
        std::copy(rValues.begin(), rValues.end(), out.As<double>());
    return out;
}

NDArray if_i64_vector(const std::vector<std::int64_t>& rValues) {
    return ra_int_array(rValues);
}

}  // namespace

ContactPairsResult contact_pairs(const Mesh& rSlave, const RegionSelector& rSlavePoints,
                                 const Mesh& rMaster, const RegionSelector& rMasterCells,
                                 const ContactPairsOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(
            "meshio++: contact_pairs: tolerance must be finite and non-negative");
    RegionSelector point_selector = rSlavePoints;
    if (point_selector.mKind < 0)
        point_selector.mKind = static_cast<std::int32_t>(RegionKind::Point);
    const std::size_t point_region_index = ra_find_region(rSlave, point_selector);
    const Region& point_region = rSlave.Region(point_region_index);
    if (point_region.mKind != RegionKind::Point)
        throw std::invalid_argument(
            "meshio++: contact_pairs: slave selector must name a Point region");
    const auto master_mask =
        ra_region_mask(rMaster, rMasterCells, RegionKind::Cell, "contact_pairs");
    const std::vector<RaFacet> master_facets = if_contact_facets(rMaster, master_mask);
    if (master_facets.empty())
        throw std::invalid_argument(
            "meshio++: contact_pairs: the master region has no queryable facets");

    double tolerance = rOptions.mTolerance;
    if (tolerance == 0.0)
        tolerance = 0.01 * if_mean_edge(rMaster, master_facets);
    const double tol_sq = tolerance * tolerance;
    std::vector<std::int64_t> slave_ids, master_cells, master_local_facets, master_subfacets,
        unmatched;
    std::vector<double> local, closest, gaps, normals;
    slave_ids.reserve(point_region.mEntries.Size());
    for (std::size_t i = 0; i < point_region.mEntries.Size(); ++i) {
        const std::int64_t node = detail::read_int(point_region.mEntries, i);
        if (node < 0 || static_cast<std::size_t>(node) >= rSlave.NumPoints())
            throw std::invalid_argument(
                "meshio++: contact_pairs: slave Point region contains an invalid point id");
        const detail::Vec3 query = ra_point(rSlave.Points(), rSlave.PointDim(), node);
        IfProjection best;
        for (const RaFacet& facet : master_facets) {
            const IfProjection hit = if_project_facet(rMaster, query, facet);
            if (!hit.mFound)
                continue;
            if (!best.mFound || hit.mDistanceSq < best.mDistanceSq ||
                (hit.mDistanceSq == best.mDistanceSq &&
                 std::tie(hit.mCell, hit.mFacet) < std::tie(best.mCell, best.mFacet)))
                best = hit;
        }
        slave_ids.push_back(node);
        if (!best.mFound || best.mDistanceSq > tol_sq) {
            master_cells.push_back(-1);
            master_local_facets.push_back(-1);
            master_subfacets.push_back(-1);
            local.insert(local.end(), {0.0, 0.0, 0.0});
            closest.insert(closest.end(), {0.0, 0.0, 0.0});
            gaps.push_back(0.0);
            normals.insert(normals.end(), {0.0, 0.0, 0.0});
            unmatched.push_back(node);
            continue;
        }
        master_cells.push_back(best.mCell);
        master_local_facets.push_back(best.mFacet);
        master_subfacets.push_back(best.mSubfacet);
        local.insert(local.end(), best.mWeights.begin(), best.mWeights.end());
        closest.insert(closest.end(), best.mPoint.begin(), best.mPoint.end());
        const double gap = detail::vec3_dot(detail::vec3_sub(query, best.mPoint), best.mNormal);
        gaps.push_back(gap);
        normals.insert(normals.end(), best.mNormal.begin(), best.mNormal.end());
    }
    if (rOptions.mRequireComplete && !unmatched.empty())
        throw std::invalid_argument("meshio++: contact_pairs: " + std::to_string(unmatched.size()) +
                                    " slave point(s) have no master facet within tolerance");

    ContactPairsResult result;
    result.mSlavePoint = if_i64_vector(slave_ids);
    result.mMasterCell = if_i64_vector(master_cells);
    result.mMasterFacet = if_i64_vector(master_local_facets);
    result.mMasterSubfacet = if_i64_vector(master_subfacets);
    result.mLocalCoordinates = if_f64_matrix(local, slave_ids.size(), 3);
    result.mClosestPoint = if_f64_matrix(closest, slave_ids.size(), 3);
    result.mGap = ra_double_array(gaps);
    result.mNormal = if_f64_matrix(normals, slave_ids.size(), 3);
    result.mUnmatched = if_i64_vector(unmatched);
    return result;
}

namespace {

std::size_t if_facet_corner_count(CellType Type, std::size_t NumNodes) {
    const int corners = detail::cell_corner_count(Type);
    return corners > 0 ? static_cast<std::size_t>(corners) : NumNodes;
}

struct SplitFacetOwner {
    RaFacet mFacet;
    std::vector<std::int64_t> mFullNodes;
    CellType mType = CellType::Custom;
};

using SplitFacetMap = std::map<std::vector<std::int64_t>, std::vector<SplitFacetOwner>>;

std::vector<std::vector<std::int64_t>> if_mesh_cell_nodes(const Mesh& rMesh) {
    std::vector<std::vector<std::int64_t>> cells;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                "meshio++: split_interface: polyhedron cells are not supported");
        for (std::size_t row = 0; row < cb.NumCells(); ++row) {
            std::vector<std::int64_t> nodes;
            if (cb.IsRagged()) {
                const std::int64_t* p = cb.Row(row);
                nodes.assign(p, p + cb.RowSize(row));
            } else {
                const std::size_t width = cb.NodesPerCell();
                const NDArray& conn = cb.Conn();
                nodes.reserve(width);
                for (std::size_t k = 0; k < width; ++k)
                    nodes.push_back(detail::read_int(conn, row * width + k));
            }
            cells.push_back(std::move(nodes));
        }
    }
    return cells;
}

SplitFacetMap if_split_facets(const Mesh& rMesh) {
    SplitFacetMap facets;
    const auto bases = detail::block_bases(rMesh);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type_name(cb.Type());
        const CellType cell_type =
            cb.IsRagged() ? CellType::Polygon : cell_type_from_name(type_name);
        if (cell_type_dimension(cell_type) != 2 && cell_type_dimension(cell_type) != 3)
            continue;
        for (std::size_t row = 0; row < cb.NumCells(); ++row) {
            const std::int64_t cell = bases[b] + static_cast<std::int64_t>(row);
            if (!cb.IsRagged() && type_name.rfind("polygon", 0) == 0) {
                const std::size_t n = cb.NodesPerCell();
                for (std::size_t edge = 0; edge < n; ++edge) {
                    const std::int64_t a = detail::read_int(cb.Conn(), row * n + edge);
                    const std::int64_t c = detail::read_int(cb.Conn(), row * n + (edge + 1) % n);
                    const std::int64_t low = std::min(a, c), high = std::max(a, c);
                    facets[{low, high}].push_back(
                        SplitFacetOwner{RaFacet{{a, c}, cell, static_cast<std::int64_t>(edge)},
                                        {a, c},
                                        CellType::Line});
                }
                continue;
            }
            const std::size_t facet_limit = cb.IsRagged() ? cb.RowSize(row) : 64;
            for (std::size_t f = 0; f < facet_limit; ++f) {
                CellType facet_type = CellType::Custom;
                std::vector<std::int64_t> full_nodes;
                if (!detail::facet_nodes(rMesh, cell, static_cast<std::int64_t>(f), facet_type,
                                         full_nodes))
                    break;
                const std::size_t corners = if_facet_corner_count(facet_type, full_nodes.size());
                if (corners < 2 || corners > full_nodes.size())
                    continue;
                std::vector<std::int64_t> corner_nodes(full_nodes.begin(),
                                                       full_nodes.begin() + corners);
                std::vector<std::int64_t> key = corner_nodes;
                std::sort(key.begin(), key.end());
                if (std::adjacent_find(key.begin(), key.end()) != key.end())
                    continue;
                facets[key].push_back(SplitFacetOwner{
                    RaFacet{std::move(corner_nodes), cell, static_cast<std::int64_t>(f)},
                    std::move(full_nodes), facet_type});
            }
        }
    }
    for (auto& [key, owners] : facets) {
        (void)key;
        std::sort(owners.begin(), owners.end(),
                  [](const SplitFacetOwner& rA, const SplitFacetOwner& rB) {
                      return std::tie(rA.mFacet.mCell, rA.mFacet.mFacet) <
                             std::tie(rB.mFacet.mCell, rB.mFacet.mFacet);
                  });
        owners.erase(std::unique(owners.begin(), owners.end(),
                                 [](const SplitFacetOwner& rA, const SplitFacetOwner& rB) {
                                     return rA.mFacet.mCell == rB.mFacet.mCell;
                                 }),
                     owners.end());
    }
    return facets;
}

struct CohesiveCell {
    std::string mType;
    std::vector<std::int64_t> mTraceA;
    std::vector<std::int64_t> mTraceB;
};

NDArray if_i64_matrix(const std::vector<std::int64_t>& rValues, std::size_t Rows, std::size_t Cols,
                      std::int64_t Fill = 0) {
    NDArray out = NDArray::Uninit(DType::Int64, {Rows, Cols});
    if (Rows == 0 || Cols == 0)
        return out;
    std::int64_t* p = out.As<std::int64_t>();
    std::fill(p, p + Rows * Cols, Fill);
    if (!rValues.empty())
        std::copy(rValues.begin(), rValues.end(), p);
    return out;
}

NDArray if_cell_block_data(const NDArray& rSource, std::size_t Rows) {
    std::vector<std::size_t> shape = rSource.Shape();
    if (shape.empty())
        shape = {Rows};
    else
        shape[0] = Rows;
    return NDArray(rSource.Dtype(), std::move(shape));
}

}  // namespace

SplitInterfaceResult split_interface(const Mesh& rMesh, const RegionSelector& rSide,
                                     const SplitInterfaceOptions& rOptions) {
    RegionSelector side_selector = rSide;
    if (side_selector.mKind < 0)
        side_selector.mKind = static_cast<std::int32_t>(RegionKind::Side);
    const std::size_t side_region_index = ra_find_region(rMesh, side_selector);
    const Region& side_region = rMesh.Region(side_region_index);
    if (side_region.mKind != RegionKind::Side)
        throw std::invalid_argument("meshio++: split_interface: selector must name a Side region");
    if (side_region.mEntries.Shape().size() != 2 || side_region.mEntries.Shape()[1] != 2)
        throw std::invalid_argument(
            "meshio++: split_interface: Side region entries must have shape (n, 2)");

    const std::vector<std::vector<std::int64_t>> cell_nodes = if_mesh_cell_nodes(rMesh);
    const SplitFacetMap facets = if_split_facets(rMesh);
    if (side_region.mEntries.Shape()[0] == 0)
        throw std::invalid_argument("meshio++: split_interface: Side region is empty");
    const auto bases = detail::block_bases(rMesh);
    std::vector<std::pair<std::int64_t, std::int64_t>> selected_entries;
    std::set<std::vector<std::int64_t>> cut_keys;
    std::set<std::int64_t> seed_cells;
    for (std::size_t i = 0; i < side_region.mEntries.Shape()[0]; ++i) {
        const std::int64_t cell = detail::read_int(side_region.mEntries, 2 * i);
        const std::int64_t facet = detail::read_int(side_region.mEntries, 2 * i + 1);
        if (cell < 0 || static_cast<std::size_t>(cell) >= cell_nodes.size())
            throw std::invalid_argument(
                "meshio++: split_interface: Side region has an invalid cell id");
        CellType facet_type = CellType::Custom;
        std::vector<std::int64_t> full_nodes;
        if (!detail::facet_nodes(rMesh, cell, facet, facet_type, full_nodes))
            throw std::invalid_argument(
                "meshio++: split_interface: Side region has an invalid local facet id");
        const std::size_t corners = if_facet_corner_count(facet_type, full_nodes.size());
        std::vector<std::int64_t> key(full_nodes.begin(), full_nodes.begin() + corners);
        std::sort(key.begin(), key.end());
        if (!facets.count(key))
            throw std::invalid_argument(
                "meshio++: split_interface: selected facet has no mesh incidence");
        selected_entries.push_back({cell, facet});
        cut_keys.insert(key);
        seed_cells.insert(cell);
    }
    std::sort(selected_entries.begin(), selected_entries.end());
    selected_entries.erase(std::unique(selected_entries.begin(), selected_entries.end()),
                           selected_entries.end());

    const std::size_t num_points = rMesh.NumPoints();
    std::vector<std::vector<std::int64_t>> incident_cells(num_points);
    for (std::size_t cell = 0; cell < cell_nodes.size(); ++cell)
        for (std::int64_t node : cell_nodes[cell]) {
            if (node < 0 || static_cast<std::size_t>(node) >= num_points)
                throw std::invalid_argument(
                    "meshio++: split_interface: input connectivity has an invalid point id");
            incident_cells[static_cast<std::size_t>(node)].push_back(
                static_cast<std::int64_t>(cell));
        }
    std::vector<std::map<std::int64_t, std::set<std::int64_t>>> node_graph(num_points);
    for (std::size_t point = 0; point < num_points; ++point)
        for (std::int64_t cell : incident_cells[point])
            node_graph[point][cell];
    for (const auto& [key, owners] : facets) {
        if (owners.size() < 2 || cut_keys.count(key))
            continue;
        for (std::size_t i = 0; i < owners.size(); ++i)
            for (std::size_t j = i + 1; j < owners.size(); ++j) {
                const auto& a = owners[i];
                const auto& b = owners[j];
                std::set<std::int64_t> nodes_b(b.mFullNodes.begin(), b.mFullNodes.end());
                std::vector<std::int64_t> common;
                for (std::int64_t node : a.mFullNodes)
                    if (nodes_b.count(node))
                        common.push_back(node);
                for (std::int64_t node : common) {
                    auto& graph = node_graph[static_cast<std::size_t>(node)];
                    graph[a.mFacet.mCell].insert(b.mFacet.mCell);
                    graph[b.mFacet.mCell].insert(a.mFacet.mCell);
                }
            }
    }

    // For every original point, label connected incident-cell fans after the
    // selected facets are removed from the local dual graph.
    std::vector<std::map<std::int64_t, std::size_t>> component_for_cell(num_points);
    std::vector<std::size_t> keeper_component(num_points, 0);
    std::vector<std::map<std::size_t, std::int64_t>> duplicate_for_component(num_points);
    std::vector<std::pair<std::int64_t, std::size_t>> duplicate_sources;
    for (std::size_t point = 0; point < num_points; ++point) {
        auto& graph = node_graph[point];
        std::set<std::int64_t> unseen;
        for (const auto& [cell, neighbours] : graph) {
            (void)neighbours;
            unseen.insert(cell);
        }
        std::vector<std::vector<std::int64_t>> components;
        while (!unseen.empty()) {
            const std::int64_t seed = *unseen.begin();
            unseen.erase(unseen.begin());
            std::vector<std::int64_t> stack{seed}, component;
            while (!stack.empty()) {
                const std::int64_t cell = stack.back();
                stack.pop_back();
                component.push_back(cell);
                for (std::int64_t neighbour : graph[cell])
                    if (unseen.erase(neighbour))
                        stack.push_back(neighbour);
            }
            std::sort(component.begin(), component.end());
            components.push_back(std::move(component));
        }
        std::sort(components.begin(), components.end(),
                  [](const auto& rA, const auto& rB) { return rA.front() < rB.front(); });
        if (components.empty())
            continue;
        std::size_t keeper = 0;
        std::size_t best_seed_count = 0;
        for (std::size_t ci = 0; ci < components.size(); ++ci) {
            std::size_t seed_count = 0;
            for (std::int64_t cell : components[ci])
                if (seed_cells.count(cell))
                    ++seed_count;
            if (seed_count > best_seed_count) {
                best_seed_count = seed_count;
                keeper = ci;
            }
        }
        keeper_component[point] = keeper;
        for (std::size_t ci = 0; ci < components.size(); ++ci)
            for (std::int64_t cell : components[ci])
                component_for_cell[point][cell] = ci;
        for (std::size_t ci = 0; ci < components.size(); ++ci)
            if (ci != keeper) {
                const std::int64_t new_id =
                    static_cast<std::int64_t>(num_points + duplicate_sources.size());
                duplicate_for_component[point][ci] = new_id;
                duplicate_sources.push_back({static_cast<std::int64_t>(point), ci});
            }
    }

    auto remap_node = [&](std::int64_t Node, std::int64_t Cell) {
        const std::size_t p = static_cast<std::size_t>(Node);
        const auto comp = component_for_cell[p].find(Cell);
        if (comp == component_for_cell[p].end())
            return Node;
        const auto copy = duplicate_for_component[p].find(comp->second);
        return copy == duplicate_for_component[p].end() ? Node : copy->second;
    };

    std::vector<CohesiveCell> cohesive;
    if (rOptions.mAddCohesive) {
        std::set<
            std::pair<std::pair<std::int64_t, std::int64_t>, std::pair<std::int64_t, std::int64_t>>>
            cohesive_interfaces;
        for (const auto& [cell, local_facet] : selected_entries) {
            CellType selected_type = CellType::Custom;
            std::vector<std::int64_t> selected_nodes;
            detail::facet_nodes(rMesh, cell, local_facet, selected_type, selected_nodes);
            std::vector<std::int64_t> key(
                selected_nodes.begin(),
                selected_nodes.begin() +
                    if_facet_corner_count(selected_type, selected_nodes.size()));
            std::sort(key.begin(), key.end());
            const auto& owners = facets.at(key);
            for (const SplitFacetOwner& other : owners) {
                if (other.mFacet.mCell == cell)
                    continue;
                auto owner_a = std::make_pair(cell, local_facet);
                auto owner_b = std::make_pair(other.mFacet.mCell, other.mFacet.mFacet);
                if (owner_b < owner_a)
                    std::swap(owner_a, owner_b);
                if (!cohesive_interfaces.insert({owner_a, owner_b}).second)
                    continue;
                const std::size_t corners =
                    if_facet_corner_count(selected_type, selected_nodes.size());
                if (corners != 2 && corners != 3 && corners != 4)
                    throw std::invalid_argument(
                        "meshio++: split_interface: cohesive cells support only line, triangle and "
                        "quad facets");
                std::vector<std::int64_t> a(selected_nodes.begin(),
                                            selected_nodes.begin() + corners);
                std::vector<std::int64_t> other_corners(
                    other.mFullNodes.begin(),
                    other.mFullNodes.begin() +
                        if_facet_corner_count(other.mType, other.mFullNodes.size()));
                if (other_corners.size() != corners)
                    throw std::invalid_argument(
                        "meshio++: split_interface: matched facet arities differ");
                std::vector<std::int64_t> b;
                b.reserve(corners);
                for (std::int64_t node : selected_nodes) {
                    const auto it = std::find(other_corners.begin(), other_corners.end(), node);
                    if (it == other_corners.end())
                        throw std::invalid_argument(
                            "meshio++: split_interface: matched facets have different corner "
                            "nodes");
                    b.push_back(node);
                }
                for (std::int64_t& node : a)
                    node = remap_node(node, cell);
                for (std::int64_t& node : b)
                    node = remap_node(node, other.mFacet.mCell);
                const std::string type = corners == 2   ? "line"
                                         : corners == 3 ? "wedge"
                                                        : "hexahedron";
                cohesive.push_back({type, std::move(a), std::move(b)});
            }
        }
        std::sort(cohesive.begin(), cohesive.end(),
                  [](const CohesiveCell& rA, const CohesiveCell& rB) {
                      return std::tie(rA.mType, rA.mTraceA, rA.mTraceB) <
                             std::tie(rB.mType, rB.mTraceA, rB.mTraceB);
                  });
        cohesive.erase(std::unique(cohesive.begin(), cohesive.end(),
                                   [](const CohesiveCell& rA, const CohesiveCell& rB) {
                                       return rA.mType == rB.mType && rA.mTraceA == rB.mTraceA &&
                                              rA.mTraceB == rB.mTraceB;
                                   }),
                       cohesive.end());
    }

    const std::vector<std::int64_t> point_sources = [&]() {
        std::vector<std::int64_t> ids(num_points);
        for (std::size_t i = 0; i < num_points; ++i)
            ids[i] = static_cast<std::int64_t>(i);
        for (const auto& [point, comp] : duplicate_sources) {
            (void)comp;
            ids.push_back(point);
        }
        return ids;
    }();
    Mesh out;
    out.AssignPoints(detail::subset_gather_rows(rMesh.Points(), point_sources));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const std::int64_t base = bases[b];
        if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows;
            rows.reserve(cb.NumCells());
            for (std::size_t row = 0; row < cb.NumCells(); ++row) {
                const std::int64_t* p = cb.Row(row);
                std::vector<std::int64_t> nodes(p, p + cb.RowSize(row));
                const std::int64_t cell = base + static_cast<std::int64_t>(row);
                for (std::int64_t& node : nodes)
                    node = remap_node(node, cell);
                rows.push_back(std::move(nodes));
            }
            out.AddPolygonBlock(type, std::move(rows));
        } else {
            const std::size_t ncell = cb.NumCells(), width = cb.NodesPerCell();
            NDArray conn = NDArray::Uninit(DType::Int64, {ncell, width});
            std::int64_t* p = conn.As<std::int64_t>();
            for (std::size_t row = 0; row < ncell; ++row)
                for (std::size_t k = 0; k < width; ++k) {
                    const std::int64_t node = detail::read_int(cb.Conn(), row * width + k);
                    p[row * width + k] = remap_node(node, base + static_cast<std::int64_t>(row));
                }
            out.AddCellBlock(type, std::move(conn));
        }
    }

    std::map<std::string, std::vector<std::size_t>> cohesive_by_type;
    for (std::size_t i = 0; i < cohesive.size(); ++i)
        cohesive_by_type[cohesive[i].mType].push_back(i);
    for (const std::string& type : {"line", "wedge", "hexahedron"}) {
        const auto it = cohesive_by_type.find(type);
        if (it == cohesive_by_type.end())
            continue;
        const std::size_t width = type == "line" ? 2 : type == "wedge" ? 6 : 8;
        NDArray conn = NDArray::Uninit(DType::Int64, {it->second.size(), width});
        std::int64_t* p = conn.As<std::int64_t>();
        for (std::size_t row = 0; row < it->second.size(); ++row) {
            const CohesiveCell& cell = cohesive[it->second[row]];
            std::copy(cell.mTraceA.begin(), cell.mTraceA.end(), p + row * width);
            if (type != "line")
                std::copy(cell.mTraceB.begin(), cell.mTraceB.end(),
                          p + row * width + cell.mTraceA.size());
        }
        out.AddCellBlock(type, std::move(conn));
    }

    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::subset_gather_rows(rMesh.PointData(name), point_sources));
    const bool has_line_cohesive = cohesive_by_type.count("line") != 0;
    for (const std::string& name : rMesh.CellDataNames()) {
        if (has_line_cohesive && name == "cohesive:trace_b")
            throw std::invalid_argument(
                "meshio++: split_interface: input already has cell_data 'cohesive:trace_b'");
        std::vector<NDArray> arrays;
        arrays.reserve(out.NumCellBlocks());
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            arrays.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
        for (const std::string& type : {"line", "wedge", "hexahedron"}) {
            const std::size_t n =
                cohesive_by_type.count(type) ? cohesive_by_type.at(type).size() : 0;
            if (n == 0)
                continue;
            const std::size_t block = arrays.size();
            const std::size_t count = static_cast<std::size_t>(out.Cells(block).NumCells());
            arrays.push_back(if_cell_block_data(rMesh.CellData(name, 0), count));
        }
        out.AddCellData(name, std::move(arrays));
    }
    if (has_line_cohesive) {
        std::vector<NDArray> arrays;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            arrays.push_back(if_i64_matrix({}, rMesh.Cells(b).NumCells(), 2, -1));
        for (const std::string& type : {"line", "wedge", "hexahedron"}) {
            const auto it = cohesive_by_type.find(type);
            if (it == cohesive_by_type.end())
                continue;
            const std::size_t rows = it->second.size();
            NDArray values = if_i64_matrix({}, rows, 2, -1);
            if (type == "line") {
                std::int64_t* p = values.As<std::int64_t>();
                for (std::size_t i = 0; i < rows; ++i) {
                    const auto& trace = cohesive[it->second[i]].mTraceB;
                    p[2 * i] = trace[0];
                    p[2 * i + 1] = trace[1];
                }
            }
            arrays.push_back(std::move(values));
        }
        out.AddCellData("cohesive:trace_b", std::move(arrays));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        Region region = rMesh.Region(i);
        if (region.mKind == RegionKind::Point) {
            std::set<std::int64_t> entries;
            for (std::size_t e = 0; e < region.mEntries.Size(); ++e)
                entries.insert(detail::read_int(region.mEntries, e));
            for (std::size_t copy = 0; copy < duplicate_sources.size(); ++copy)
                if (entries.count(duplicate_sources[copy].first))
                    entries.insert(static_cast<std::int64_t>(num_points + copy));
            std::vector<std::int64_t> flat(entries.begin(), entries.end());
            region.mEntries = ra_int_array(flat);
        }
        out.AddRegion(std::move(region));
    }
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        out.AddPropertySet(rMesh.GetPropertySet(i));

    SplitInterfaceResult result;
    result.mMesh = std::move(out);
    result.mNumDuplicatedPoints = static_cast<std::int64_t>(duplicate_sources.size());
    result.mNumCohesiveCells = static_cast<std::int64_t>(cohesive.size());
    return result;
}

}  // namespace meshioplusplus
