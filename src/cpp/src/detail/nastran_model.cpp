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
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Nastran numbers the hex20/wedge15 mid-side nodes bottom, vertical, top;
// meshio++ (VTK) numbers them bottom, top, vertical (the bulk reader's tables).
constexpr int kNmHexa20[20] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                               10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int kNmPenta15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};

constexpr NastranCardSpec kNmCards[] = {
    {"CBAR", "line", 2, nullptr, 0, nullptr},
    {"CBEAM", "line", 2, nullptr, 0, nullptr},
    {"CBUSH", "line", 2, nullptr, 0, nullptr},
    // Springs and dampers between two points, or grounded (the second is 0).
    {"CDAMP1", "vertex", 1, "line", 2, nullptr},
    {"CDAMP2", "vertex", 1, "line", 2, nullptr},
    {"CELAS1", "vertex", 1, "line", 2, nullptr},
    {"CELAS2", "vertex", 1, "line", 2, nullptr},
    {"CHEXA", "hexahedron", 8, "hexahedron20", 20, kNmHexa20},
    {"CONM2", "vertex", 1, nullptr, 0, nullptr},
    {"CONROD", "line", 2, nullptr, 0, nullptr},
    {"CPENTA", "wedge", 6, "wedge15", 15, kNmPenta15},
    {"CPYRAM", "pyramid", 5, "pyramid13", 13, nullptr},
    {"CQUAD", "quad", 4, "quad9", 9, nullptr},
    {"CQUAD4", "quad", 4, nullptr, 0, nullptr},
    {"CQUAD8", "quad", 4, "quad8", 8, nullptr},
    {"CQUADR", "quad", 4, nullptr, 0, nullptr},
    {"CROD", "line", 2, nullptr, 0, nullptr},
    {"CSHEAR", "quad", 4, nullptr, 0, nullptr},
    {"CTETRA", "tetra", 4, "tetra10", 10, nullptr},
    {"CTRIA3", "triangle", 3, nullptr, 0, nullptr},
    {"CTRIA6", "triangle", 3, "triangle6", 6, nullptr},
    {"CTRIAR", "triangle", 3, nullptr, 0, nullptr},
    {"CTUBE", "line", 2, nullptr, 0, nullptr},
    {"CVISC", "line", 2, nullptr, 0, nullptr},
    {"PLOTEL", "line", 2, nullptr, 0, nullptr},
};

NDArray nm_int_array(const std::vector<std::int64_t>& rV) {
    NDArray a(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

struct NmBlock {
    std::string mCard;
    std::string mType;
    std::size_t mNodes;
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mEid;
    std::vector<std::int64_t> mPid;
};

}  // namespace

const NastranCardSpec* nastran_card_spec(std::string_view Card) {
    for (const NastranCardSpec& spec : kNmCards)
        if (Card == spec.mCard)
            return &spec;
    return nullptr;
}

NastranCells nastran_add_cells(Mesh& rMesh, const std::vector<NastranCardRows>& rCards,
                               const std::unordered_map<std::int64_t, std::size_t>& rGridIndex,
                               const std::unordered_set<std::int64_t>& rScalarPoints,
                               const std::map<std::int64_t, std::string>& rPtype,
                               const std::string& rWho) {
    std::vector<NmBlock> blocks;
    std::size_t dropped = 0;
    for (const NastranCardRows& card : rCards) {
        const NastranCardSpec* spec = nastran_card_spec(card.mCard);
        if (spec == nullptr)
            continue;
        const std::size_t width = card.mWidth;
        const std::size_t n = card.mEid.size();
        NmBlock linear{card.mCard, spec->mLinear, spec->mLinearNodes, {}, {}, {}};
        NmBlock quadratic{
            card.mCard, spec->mQuadratic ? spec->mQuadratic : "", spec->mQuadraticNodes, {}, {},
            {}};
        std::size_t partial = 0;
        // A spring or damper grounded at its first end: the other end leads.
        const bool grounded_first =
            spec->mLinearNodes == 1 && spec->mQuadraticNodes == 2 && width >= 2;
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t* row = card.mNodes.data() + i * width;
            const std::int64_t swapped[2] = {grounded_first ? row[1] : 0, 0};
            if (grounded_first && row[0] == 0)
                row = swapped;
            bool quad = false;
            if (spec->mQuadratic != nullptr && width >= spec->mQuadraticNodes) {
                std::size_t given = 0;
                for (std::size_t k = spec->mLinearNodes; k < spec->mQuadraticNodes; ++k)
                    given += row[k] != 0 ? 1 : 0;
                quad = given == spec->mQuadraticNodes - spec->mLinearNodes;
                partial += (given != 0 && !quad) ? 1 : 0;
            }
            NmBlock& b = quad ? quadratic : linear;
            std::vector<std::int64_t> conn(b.mNodes);
            bool ok = true;
            for (std::size_t k = 0; k < b.mNodes && ok; ++k) {
                const std::size_t src = (quad && spec->mPermutation)
                                            ? static_cast<std::size_t>(spec->mPermutation[k])
                                            : k;
                const auto it = rGridIndex.find(row[src]);
                if (it == rGridIndex.end()) {
                    if (rScalarPoints.count(row[src]) == 0)
                        throw ReadError(rWho + ": " + card.mCard + " " +
                                        std::to_string(card.mEid[i]) + " references " +
                                        (row[src] == 0
                                             ? std::string("no node")
                                             : "undefined GRID " + std::to_string(row[src])) +
                                        " as its node " + std::to_string(src + 1));
                    ok = false;  // a scalar point, not a GRID
                } else {
                    conn[k] = static_cast<std::int64_t>(it->second);
                }
            }
            if (!ok) {
                ++dropped;
                continue;
            }
            b.mConn.insert(b.mConn.end(), conn.begin(), conn.end());
            b.mEid.push_back(card.mEid[i]);
            b.mPid.push_back(card.mPid[i]);
        }
        if (partial != 0)
            log::warn("{}: {} {} element(s) have only some mid-side nodes; read as {}", rWho,
                      partial, card.mCard, spec->mLinear);
        if (!linear.mEid.empty())
            blocks.push_back(std::move(linear));
        if (!quadratic.mEid.empty())
            blocks.push_back(std::move(quadratic));
    }
    if (dropped != 0)
        log::warn("{}: skipped {} element(s) that connect scalar points", rWho, dropped);

    NastranCells out;
    {
        std::vector<NDArray> eids;
        std::vector<NDArray> pids;
        std::size_t shared = 0;
        for (NmBlock& b : blocks) {
            out.mOffsets.push_back(out.mNumCells);
            out.mSizes.push_back(b.mEid.size());
            if (b.mCard != "CONM2")
                for (std::size_t i = 0; i < b.mEid.size(); ++i)
                    shared += out.mCellIndex.emplace(b.mEid[i], out.mNumCells + i).second ? 0 : 1;
            out.mNumCells += b.mEid.size();
            NDArray conn(DType::Int64, {b.mEid.size(), b.mNodes});
            std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
            rMesh.AddCellBlock(b.mType, std::move(conn));
            eids.push_back(nm_int_array(b.mEid));
            pids.push_back(nm_int_array(b.mPid));
        }
        if (shared != 0)
            log::warn(
                "{}: {} element id(s) are used by more than one card; their element results "
                "go to the first",
                rWho, shared);
        if (!blocks.empty()) {
            rMesh.AddCellData("nastran:eid", std::move(eids));
            rMesh.AddCellData("nastran:pid", std::move(pids));
        }
    }

    std::map<std::int64_t, std::pair<std::vector<std::int64_t>, int>> by_pid;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        const int dim = cell_type_dimension(cell_type_from_name(blocks[b].mType));
        for (std::size_t i = 0; i < blocks[b].mPid.size(); ++i) {
            const std::int64_t p = blocks[b].mPid[i];
            if (p <= 0)
                continue;
            auto& entry = by_pid[p];
            if (entry.first.empty())
                entry.second = dim;
            entry.first.push_back(static_cast<std::int64_t>(out.mOffsets[b] + i));
            entry.second = std::max(entry.second, dim);
        }
    }
    for (auto& [p, entry] : by_pid) {
        const auto it = rPtype.find(p);
        const std::string name =
            (it != rPtype.end() ? it->second : std::string("PID")) + "_" + std::to_string(p);
        rMesh.AddRegion(Region(name, RegionKind::Cell, entry.second, p, nm_int_array(entry.first)));
    }
    return out;
}

namespace {

// Degrees to radians; one constant so both engines round the same way.
constexpr double kNmDegree = 3.14159265358979323846 / 180.0;

void nm_local_to_cartesian(int Type, const double* pLocal, double* pOut) {
    if (Type == 2) {
        const double t = pLocal[1] * kNmDegree;
        pOut[0] = pLocal[0] * std::cos(t);
        pOut[1] = pLocal[0] * std::sin(t);
        pOut[2] = pLocal[2];
    } else if (Type == 3) {
        const double t = pLocal[1] * kNmDegree;
        const double f = pLocal[2] * kNmDegree;
        pOut[0] = pLocal[0] * std::sin(t) * std::cos(f);
        pOut[1] = pLocal[0] * std::sin(t) * std::sin(f);
        pOut[2] = pLocal[0] * std::cos(t);
    } else {
        pOut[0] = pLocal[0];
        pOut[1] = pLocal[1];
        pOut[2] = pLocal[2];
    }
}

// A system from its origin A, a point B on +z and a point C in the xz plane,
// all in basic; false when they do not span one.
bool nm_system_from_points(int Type, const double* pA, const double* pB, const double* pC,
                           NastranCoordSystems::System& rOut) {
    const double z[3] = {pB[0] - pA[0], pB[1] - pA[1], pB[2] - pA[2]};
    const double nz = std::sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
    if (!(nz > 0.0))
        return false;
    const double ez[3] = {z[0] / nz, z[1] / nz, z[2] / nz};
    const double v[3] = {pC[0] - pA[0], pC[1] - pA[1], pC[2] - pA[2]};
    const double y[3] = {ez[1] * v[2] - ez[2] * v[1], ez[2] * v[0] - ez[0] * v[2],
                         ez[0] * v[1] - ez[1] * v[0]};
    const double ny = std::sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    if (!(ny > 0.0))
        return false;
    const double ey[3] = {y[0] / ny, y[1] / ny, y[2] / ny};
    const double ex[3] = {ey[1] * ez[2] - ey[2] * ez[1], ey[2] * ez[0] - ey[0] * ez[2],
                          ey[0] * ez[1] - ey[1] * ez[0]};
    rOut.mType = Type;
    for (int k = 0; k < 3; ++k) {
        rOut.mOrigin[k] = pA[k];
        rOut.mAxes[k] = ex[k];
        rOut.mAxes[3 + k] = ey[k];
        rOut.mAxes[6 + k] = ez[k];
    }
    return true;
}

}  // namespace

NastranCoordSystems::NastranCoordSystems() {
    mSystems[0] = System{};
}

void NastranCoordSystems::ToBasic(std::int64_t Cid, const double* pLocal, double* pBasic) const {
    const System& s = mSystems.at(Cid);
    double c[3];
    nm_local_to_cartesian(s.mType, pLocal, c);
    for (int k = 0; k < 3; ++k)
        pBasic[k] =
            s.mOrigin[k] + s.mAxes[k] * c[0] + s.mAxes[3 + k] * c[1] + s.mAxes[6 + k] * c[2];
}

void NastranCoordSystems::VectorToBasic(std::int64_t Cid, const double* pBasicPoint,
                                        double* pVector) const {
    const System& s = mSystems.at(Cid);
    double w[3] = {pVector[0], pVector[1], pVector[2]};
    if (s.mType == 2 || s.mType == 3) {
        const double d[3] = {pBasicPoint[0] - s.mOrigin[0], pBasicPoint[1] - s.mOrigin[1],
                             pBasicPoint[2] - s.mOrigin[2]};
        double l[3];
        for (int j = 0; j < 3; ++j)
            l[j] = s.mAxes[3 * j] * d[0] + s.mAxes[3 * j + 1] * d[1] + s.mAxes[3 * j + 2] * d[2];
        const double ph = std::atan2(l[1], l[0]);
        const double cp = std::cos(ph);
        const double sp = std::sin(ph);
        if (s.mType == 2) {
            w[0] = pVector[0] * cp - pVector[1] * sp;
            w[1] = pVector[0] * sp + pVector[1] * cp;
        } else {
            const double th = std::atan2(std::sqrt(l[0] * l[0] + l[1] * l[1]), l[2]);
            const double ct = std::cos(th);
            const double st = std::sin(th);
            const double er[3] = {st * cp, st * sp, ct};
            const double et[3] = {ct * cp, ct * sp, -st};
            const double ep[3] = {-sp, cp, 0.0};
            for (int k = 0; k < 3; ++k)
                w[k] = pVector[0] * er[k] + pVector[1] * et[k] + pVector[2] * ep[k];
        }
    }
    for (int k = 0; k < 3; ++k)
        pVector[k] = s.mAxes[k] * w[0] + s.mAxes[3 + k] * w[1] + s.mAxes[6 + k] * w[2];
}

NastranCoordSystems nastran_apply_frames(Mesh& rMesh, const std::vector<NastranCoordCard>& rCards,
                                         const std::vector<std::int64_t>& rIds,
                                         const std::vector<std::int64_t>& rCp,
                                         const std::vector<std::int64_t>& rCd,
                                         const std::string& rWho) {
    NastranCoordSystems systems;
    const std::size_t n = rIds.size();
    const bool any_cp = std::any_of(rCp.begin(), rCp.end(), [](std::int64_t c) { return c != 0; });
    const bool any_cd = std::any_of(rCd.begin(), rCd.end(), [](std::int64_t c) { return c != 0; });
    if (!any_cp && !any_cd)
        return systems;

    // First definition of each CID wins; later ones are reported.
    std::map<std::int64_t, const NastranCoordCard*> pending;
    std::size_t duplicates = 0;
    for (const NastranCoordCard& c : rCards) {
        if (c.mCid <= 0)
            continue;
        if (!pending.emplace(c.mCid, &c).second)
            ++duplicates;
    }
    if (duplicates != 0)
        log::warn(
            "{}: {} coordinate system(s) are defined more than once; the first definition "
            "is used",
            rWho, duplicates);

    NDArray points = rMesh.Points();  // deep copy: moved to basic below, then reassigned
    double* xyz = points.As<double>();
    std::unordered_map<std::int64_t, std::size_t> index;
    index.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        index.emplace(rIds[i], i);
    std::vector<char> resolved(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        resolved[i] = rCp[i] == 0 ? 1 : 0;

    // A CORD2 needs its reference system, a CORD1 its three GRIDs, and a GRID
    // its CP system: resolve whatever is ready until nothing changes.
    std::set<std::int64_t> degenerate;
    for (bool progress = true; progress;) {
        progress = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (resolved[i] || !systems.Has(rCp[i]))
                continue;
            double local[3] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
            systems.ToBasic(rCp[i], local, xyz + 3 * i);
            resolved[i] = 1;
            progress = true;
        }
        for (auto it = pending.begin(); it != pending.end();) {
            const NastranCoordCard& c = *it->second;
            double a[3], b[3], p[3];
            bool ready = false;
            if (c.mByGrids) {
                std::size_t g[3];
                ready = true;
                for (int k = 0; k < 3 && ready; ++k) {
                    const auto found = index.find(c.mGrids[k]);
                    ready = found != index.end() && resolved[found->second];
                    if (ready)
                        g[k] = found->second;
                }
                if (ready)
                    for (int k = 0; k < 3; ++k) {
                        a[k] = xyz[3 * g[0] + k];
                        b[k] = xyz[3 * g[1] + k];
                        p[k] = xyz[3 * g[2] + k];
                    }
            } else if (systems.Has(c.mRid)) {
                ready = true;
                systems.ToBasic(c.mRid, c.mAbc, a);
                systems.ToBasic(c.mRid, c.mAbc + 3, b);
                systems.ToBasic(c.mRid, c.mAbc + 6, p);
            }
            if (!ready) {
                ++it;
                continue;
            }
            NastranCoordSystems::System sys;
            if (nm_system_from_points(c.mType, a, b, p, sys))
                systems.Add(c.mCid, sys);
            else
                degenerate.insert(c.mCid);
            it = pending.erase(it);
            progress = true;
        }
    }
    if (!degenerate.empty())
        log::warn(
            "{}: {} coordinate system(s) have coincident or collinear defining points and "
            "are ignored",
            rWho, degenerate.size());

    std::size_t moved = 0, kept = 0;
    std::set<std::int64_t> missing;
    for (std::size_t i = 0; i < n; ++i) {
        if (rCp[i] == 0)
            continue;
        if (resolved[i]) {
            ++moved;
        } else {
            ++kept;
            missing.insert(rCp[i]);
        }
    }
    if (kept != 0) {
        std::string ids;
        for (std::int64_t c : missing)
            ids += (ids.empty() ? "" : ", ") + std::to_string(c);
        log::warn(
            "{}: {} GRID(s) are in coordinate system(s) {} that cannot be resolved; their "
            "coordinates are kept as written",
            rWho, kept, ids);
    }
    if (moved != 0) {
        log::info("{}: moved {} GRID(s) from local coordinate systems to basic", rWho, moved);
        rMesh.AssignPoints(std::move(points));
    }
    if (any_cp)
        rMesh.AddPointData("nastran:cp", nm_int_array(rCp));
    if (any_cd) {
        std::size_t unresolved = 0;
        for (std::int64_t c : rCd)
            unresolved += (c > 0 && !systems.Has(c)) ? 1 : 0;
        if (unresolved != 0)
            log::warn(
                "{}: {} GRID(s) have an output system (CD) that cannot be resolved; their "
                "results stay in it",
                rWho, unresolved);
        rMesh.AddPointData("nastran:cd", nm_int_array(rCd));
    }
    return systems;
}

bool nastran_rotate_to_basic(const NastranCoordSystems& rSystems, std::int64_t Cd,
                             const double* pBasicPoint, double* pValues, std::size_t Count) {
    if (Cd <= 0 || !rSystems.Has(Cd))
        return false;
    for (std::size_t k = 0; k < Count; ++k)
        rSystems.VectorToBasic(Cd, pBasicPoint, pValues + 3 * k);
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus
