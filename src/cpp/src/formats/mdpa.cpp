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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/backends/kratos_names.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Lexing helpers. Every anonymous-namespace symbol here is `mdpa_`-prefixed
// because the single-header amalgamation concatenates all of src/cpp/src.
// ---------------------------------------------------------------------------

std::string mdpa_strip(const std::string& rS) {
    const std::size_t b = rS.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    const std::size_t e = rS.find_last_not_of(" \t\r\n");
    return rS.substr(b, e - b + 1);
}

/// The line with any `//` comment removed, then stripped.
std::string mdpa_clean(const std::string& rS) {
    const std::size_t c = rS.find("//");
    return mdpa_strip(c == std::string::npos ? rS : rS.substr(0, c));
}

std::vector<std::string> mdpa_tokens(const std::string& rS) {
    std::vector<std::string> out;
    std::istringstream iss(rS);
    std::string t;
    while (iss >> t)
        out.push_back(t);
    return out;
}

bool mdpa_starts_with(const std::string& rS, const std::string& rPrefix) {
    return rS.size() >= rPrefix.size() && rS.compare(0, rPrefix.size(), rPrefix) == 0;
}

bool mdpa_parse_int(const std::string& rS, std::int64_t& rOut) {
    if (rS.empty())
        return false;
    char* end = nullptr;
    const long long v = std::strtoll(rS.c_str(), &end, 10);
    if (end != rS.c_str() + rS.size())
        return false;
    rOut = static_cast<std::int64_t>(v);
    return true;
}

bool mdpa_parse_double(const std::string& rS, double& rOut) {
    if (rS.empty())
        return false;
    char* end = nullptr;
    const double v = std::strtod(rS.c_str(), &end);
    if (end != rS.c_str() + rS.size())
        return false;
    rOut = v;
    return true;
}

// ---------------------------------------------------------------------------
// Kratos <-> VTK node orderings.
//
// `mdpa_kratos_node_order(type)[i]` is the meshio/VTK slot the node in Kratos
// slot `i` belongs in. Reading therefore scatters (`out[table[i]] = in[i]`) and
// writing gathers (`out[j] = in[table[j]]`) — the exact inverse pair the Python
// reference builds with `np.argsort`.
// ---------------------------------------------------------------------------

const std::vector<int>& mdpa_kratos_node_order(CellType type) {
    static const std::vector<int> h20 = {0,  1, 2,  3,  4,  5,  6,  7,  8,  11,
                                         10, 9, 16, 19, 18, 17, 12, 13, 14, 15};
    static const std::vector<int> h27 = {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 10, 9,  16, 19,
                                         18, 17, 12, 15, 14, 13, 20, 23, 21, 24, 22, 25, 26};
    static const std::vector<int> none;
    if (type == CellType::Hexahedron20)
        return h20;
    if (type == CellType::Hexahedron27)
        return h27;
    return none;
}

/**
 * @brief Resolve a Kratos entity name to a meshio cell type.
 *
 * Exact lookup first (`kratos_names.hpp` owns the tables), then the longest
 * suffix that resolves — which is what makes application-specific names such
 * as `SmallDisplacementElement3D4N` work, mirroring the Python reference's
 * substring fallback for the cases that occur in practice.
 */
CellType mdpa_entity_cell_type(const std::string& rName) {
    if (rName.empty())
        return CellType::Custom;
    CellType t = cell_type_from_kratos_name(rName);
    if (t != CellType::Custom)
        return t;
    for (std::size_t i = 1; i + 1 < rName.size(); ++i) {
        t = cell_type_from_kratos_name(rName.substr(i));
        if (t != CellType::Custom)
            return t;
    }
    return CellType::Custom;
}

// ---------------------------------------------------------------------------
// Reader staging
// ---------------------------------------------------------------------------

struct MdpaBlock {
    std::string mType;
    std::size_t mNodes = 0;
    bool mIsCondition = false;
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mProps;
    std::size_t mCount = 0;
};

/// One parsed row of a `NodalData` / `ElementalData` / `ConditionalData` block.
struct MdpaDataRow {
    std::size_t mBlock = 0;  ///< cell block index (0 for nodal data)
    std::size_t mRow = 0;    ///< point index, or row within the cell block
    int mFixed = -1;         ///< the optional leading 0/1 column, -1 if absent
    std::vector<double> mValues;
};

/// A cursor over the file's lines, so every block parser advances one index.
struct MdpaCursor {
    const std::vector<std::string>* mpLines = nullptr;
    std::size_t mIndex = 0;

    bool Done() const { return mIndex >= mpLines->size(); }
    const std::string& Next() { return (*mpLines)[mIndex++]; }
};

/// Consume the rest of a block, throwing if it carries anything but comments.
void mdpa_expect_empty_block(MdpaCursor& rCur, const std::string& rEnd, const std::string& rWhat) {
    while (!rCur.Done()) {
        const std::string line = mdpa_clean(rCur.Next());
        if (line.empty())
            continue;
        if (line == rEnd)
            return;
        throw ReadError("MDPA: " + rWhat +
                        " is not supported by the C++ reader (offending line: '" + line + "')");
    }
    throw ReadError("MDPA: EOF before '" + rEnd + "'");
}

/**
 * @brief Parse a `*Data` block body, mirroring `_parse_generic_data_block`.
 *
 * @param rCur cursor positioned just after the `Begin ...` header
 * @param rEnd the terminating token, e.g. `"End NodalData"`
 * @param nodal whether the optional leading `fixed` column may appear
 * @param rResolve maps a 1-based entity id to `(block, row)`; returns false to
 *        skip the line (unknown id)
 * @param rRows out: the parsed rows
 * @param rHasFixed out: whether any row carried a `fixed` column
 * @return the number of value components (0 for a membership flag), or -1 when
 *         the block held no usable line
 */
int mdpa_parse_data_block(
    MdpaCursor& rCur, const std::string& rEnd, bool nodal,
    const std::function<bool(std::int64_t, std::size_t&, std::size_t&)>& rResolve,
    std::vector<MdpaDataRow>& rRows, bool& rHasFixed) {
    int nc = -1;
    bool terminated = false;
    while (!rCur.Done()) {
        const std::string raw = mdpa_strip(rCur.Next());
        if (raw == rEnd) {
            terminated = true;
            break;
        }
        const std::string line = mdpa_clean(raw);
        if (line.empty())
            continue;
        if (line == rEnd) {
            terminated = true;
            break;
        }
        const std::vector<std::string> toks = mdpa_tokens(line);
        std::int64_t id = 0;
        if (!mdpa_parse_int(toks[0], id)) {
            log::warn("mdpa: skipping data line with non-integer id: {}", line);
            continue;
        }
        std::size_t block = 0, row = 0;
        if (!rResolve(id, block, row)) {
            log::warn("mdpa: skipping data for unknown entity id {}", id);
            continue;
        }

        // Split the remaining tokens into an optional `fixed` flag and values.
        const std::size_t rest = toks.size() - 1;
        int fixed = -1;
        std::size_t first_value = 1;
        std::int64_t flag = 0;
        const bool flag_like =
            nodal && rest > 0 && mdpa_parse_int(toks[1], flag) && (flag == 0 || flag == 1);
        if (nc < 0) {  // first usable line defines the layout
            if (flag_like && rest > 1) {
                fixed = static_cast<int>(flag);
                first_value = 2;
            }
            nc = static_cast<int>(toks.size() - first_value);
        } else {
            const std::size_t want = static_cast<std::size_t>(nc);
            if (flag_like && rest == want + 1) {
                fixed = static_cast<int>(flag);
                first_value = 2;
            } else if (rest != want) {
                log::warn("mdpa: skipping data line with {} values (expected {}): {}", rest, want,
                          line);
                continue;
            }
        }
        if (fixed >= 0)
            rHasFixed = true;

        MdpaDataRow r;
        r.mBlock = block;
        r.mRow = row;
        r.mFixed = fixed;
        bool ok = true;
        for (std::size_t j = first_value; j < toks.size(); ++j) {
            double v = 0.0;
            if (!mdpa_parse_double(toks[j], v)) {
                ok = false;
                break;
            }
            r.mValues.push_back(v);
        }
        if (!ok) {
            log::warn("mdpa: skipping data line with non-numeric value: {}", line);
            continue;
        }
        rRows.push_back(std::move(r));
    }
    if (!terminated)
        throw ReadError("MDPA: EOF before '" + rEnd + "'");
    return nc;
}

/// Build one data array of @p count rows and @p nc components from @p rRows.
NDArray mdpa_data_array(const std::vector<MdpaDataRow>& rRows, std::size_t block, std::size_t count,
                        int nc) {
    if (nc == 0) {  // membership flag
        NDArray a(DType::Int64, {count});
        std::int64_t* p = a.As<std::int64_t>();
        for (std::size_t i = 0; i < count; ++i)
            p[i] = 0;
        for (const auto& r : rRows)
            if (r.mBlock == block && r.mRow < count)
                p[r.mRow] = 1;
        return a;
    }
    const std::size_t k = static_cast<std::size_t>(nc);
    std::vector<std::size_t> shape;
    if (k == 1)
        shape = {count};
    else
        shape = {count, k};
    NDArray a(DType::Float64, shape);
    double* p = a.As<double>();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < count * k; ++i)
        p[i] = nan;
    for (const auto& r : rRows) {
        if (r.mBlock != block || r.mRow >= count || r.mValues.size() != k)
            continue;
        for (std::size_t j = 0; j < k; ++j)
            p[r.mRow * k + j] = r.mValues[j];
    }
    return a;
}

}  // namespace

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

Mesh read_mdpa(const std::string& rPath) {
    std::ifstream in(rPath);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    std::vector<std::string> lines;
    std::string l;
    while (std::getline(in, l))
        lines.push_back(l);

    MdpaCursor cur{&lines, 0};

    std::vector<double> coords;  // flat (n, 3)
    std::size_t num_points = 0;
    std::vector<MdpaBlock> blocks;
    std::unordered_map<std::int64_t, std::pair<std::size_t, std::size_t>> element_ids,
        condition_ids;
    std::map<std::string, NDArray> field_data;

    // Staged data arrays, materialized after every block is known.
    struct StagedData {
        std::string mName;
        int mComponents = 0;
        bool mHasFixed = false;
        std::vector<MdpaDataRow> mRows;
    };
    std::vector<StagedData> point_data, cell_data;

    // SubModelParts, accumulated by (hierarchical) name.
    struct StagedSmp {
        std::vector<std::int64_t> mNodes;
        std::vector<std::pair<std::size_t, std::size_t>> mCells;
    };
    std::map<std::string, StagedSmp> smps;
    std::vector<std::string> smp_stack;

    auto smp_name = [&]() -> std::string {
        std::string out;
        for (std::size_t i = 0; i < smp_stack.size(); ++i) {
            if (i)
                out += "/";
            out += smp_stack[i];
        }
        return out;
    };

    auto read_id_list = [&](const std::string& rEnd, std::vector<std::int64_t>& rOut) {
        while (!cur.Done()) {
            const std::string line = mdpa_clean(cur.Next());
            if (line.empty())
                continue;
            if (line == rEnd)
                return;
            std::int64_t v = 0;
            if (!mdpa_parse_int(line, v)) {
                log::warn("mdpa: skipping non-integer id in {}: {}", rEnd, line);
                continue;
            }
            rOut.push_back(v);
        }
        throw ReadError("MDPA: EOF before '" + rEnd + "'");
    };

    while (!cur.Done()) {
        const std::string raw = mdpa_strip(cur.Next());
        const std::string line = mdpa_clean(raw);
        if (line.empty())
            continue;

        if (line == "Begin ModelPartData") {
            while (true) {
                if (cur.Done())
                    throw ReadError("MDPA: EOF before 'End ModelPartData'");
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End ModelPartData")
                    break;
                const std::vector<std::string> t = mdpa_tokens(e);
                if (t.size() < 2) {
                    log::warn("mdpa: skipping malformed ModelPartData line: {}", e);
                    continue;
                }
                double v = 0.0;
                if (t.size() != 2 || !mdpa_parse_double(t[1], v))
                    throw ReadError("MDPA: non-numeric ModelPartData value for '" + t[0] +
                                    "' is not supported by the C++ reader");
                NDArray a(DType::Float64, {1});
                a.As<double>()[0] = v;
                field_data[t[0]] = std::move(a);
            }
        } else if (line == "Begin Nodes") {
            if (num_points)
                throw ReadError("MDPA: more than one Nodes block");
            bool terminated = false;
            while (!cur.Done()) {
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == "End Nodes") {
                    terminated = true;
                    break;
                }
                const std::vector<std::string> t = mdpa_tokens(e);
                if (t.size() < 3)
                    throw ReadError("MDPA: node line with fewer than 3 coordinates: " + e);
                if (t.size() >= 4) {
                    std::int64_t id = 0;
                    if (!mdpa_parse_int(t[0], id))
                        throw ReadError("MDPA: non-integer node id: " + e);
                    if (id != static_cast<std::int64_t>(num_points) + 1)
                        throw ReadError(
                            "MDPA: non-sequential node ids are not supported by the C++ reader "
                            "(expected " +
                            std::to_string(num_points + 1) + ", got " + std::to_string(id) + ")");
                }
                for (std::size_t c = t.size() - 3; c < t.size(); ++c) {
                    double v = 0.0;
                    if (!mdpa_parse_double(t[c], v))
                        throw ReadError("MDPA: non-numeric node coordinate: " + e);
                    coords.push_back(v);
                }
                ++num_points;
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before 'End Nodes'");
        } else if (mdpa_starts_with(line, "Begin Elements") ||
                   mdpa_starts_with(line, "Begin Conditions")) {
            const bool is_condition = mdpa_starts_with(line, "Begin Conditions");
            const std::string end_token = is_condition ? "End Conditions" : "End Elements";
            const std::vector<std::string> head = mdpa_tokens(line);
            const std::string entity_name = head.size() >= 3 ? head[2] : std::string();
            const CellType type = mdpa_entity_cell_type(entity_name);
            if (type == CellType::Custom)
                throw ReadError("MDPA: unknown Kratos entity name '" + entity_name + "'");
            const int nn = cell_type_num_nodes(type);
            if (nn <= 0)
                throw ReadError("MDPA: entity '" + entity_name +
                                "' maps to a variable-node-count cell type");
            const std::string& type_name = cell_type_name(type);
            const std::vector<int>& order = mdpa_kratos_node_order(type);

            bool terminated = false;
            while (!cur.Done()) {
                const std::string e = mdpa_clean(cur.Next());
                if (e.empty())
                    continue;
                if (e == end_token) {
                    terminated = true;
                    break;
                }
                if (mdpa_starts_with(e, "End "))
                    throw ReadError("MDPA: expected '" + end_token + "', got '" + e + "'");
                const std::vector<std::string> t = mdpa_tokens(e);
                if (static_cast<int>(t.size()) != nn + 2)
                    throw ReadError("MDPA: " + entity_name + " row with " +
                                    std::to_string(t.size() >= 2 ? t.size() - 2 : 0) +
                                    " nodes (expected " + std::to_string(nn) + "): " + e);
                std::int64_t id = 0, prop = 0;
                if (!mdpa_parse_int(t[0], id) || !mdpa_parse_int(t[1], prop))
                    throw ReadError("MDPA: non-integer id/property in: " + e);

                if (blocks.empty() || blocks.back().mType != type_name ||
                    blocks.back().mIsCondition != is_condition) {
                    MdpaBlock b;
                    b.mType = type_name;
                    b.mNodes = static_cast<std::size_t>(nn);
                    b.mIsCondition = is_condition;
                    blocks.push_back(std::move(b));
                }
                MdpaBlock& blk = blocks.back();
                const std::size_t base = blk.mConn.size();
                blk.mConn.resize(base + static_cast<std::size_t>(nn));
                for (int j = 0; j < nn; ++j) {
                    std::int64_t node = 0;
                    if (!mdpa_parse_int(t[static_cast<std::size_t>(j) + 2], node))
                        throw ReadError("MDPA: non-integer node id in: " + e);
                    const std::size_t slot =
                        order.empty()
                            ? static_cast<std::size_t>(j)
                            : static_cast<std::size_t>(order[static_cast<std::size_t>(j)]);
                    blk.mConn[base + slot] = node - 1;
                }
                blk.mProps.push_back(prop);
                const std::size_t row = blk.mCount++;
                auto& id_map = is_condition ? condition_ids : element_ids;
                id_map[id] = {blocks.size() - 1, row};
            }
            if (!terminated)
                throw ReadError("MDPA: EOF before '" + end_token + "'");
        } else if (mdpa_starts_with(line, "Begin Properties")) {
            mdpa_expect_empty_block(cur, "End Properties", "a non-empty Properties block");
        } else if (mdpa_starts_with(line, "Begin NodalData")) {
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed NodalData header: " + line);
            if (num_points == 0)
                throw ReadError("MDPA: NodalData before Nodes");
            std::string name = head[2];
            const std::size_t br = name.find('[');
            if (br != std::string::npos)
                name = name.substr(0, br);
            StagedData sd;
            sd.mName = name;
            const int nc = mdpa_parse_data_block(
                cur, "End NodalData", /*nodal=*/true,
                [&](std::int64_t id, std::size_t& block, std::size_t& row) {
                    if (id < 1 || static_cast<std::size_t>(id) > num_points)
                        return false;
                    block = 0;
                    row = static_cast<std::size_t>(id) - 1;
                    return true;
                },
                sd.mRows, sd.mHasFixed);
            if (nc < 0) {
                log::warn("mdpa: NodalData block '{}' held no usable line", name);
                continue;
            }
            sd.mComponents = nc;
            point_data.push_back(std::move(sd));
        } else if (mdpa_starts_with(line, "Begin ElementalData") ||
                   mdpa_starts_with(line, "Begin ConditionalData")) {
            const bool elemental = mdpa_starts_with(line, "Begin ElementalData");
            const std::string end_token = elemental ? "End ElementalData" : "End ConditionalData";
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed " + end_token.substr(4) + " header: " + line);
            std::string name = head[2];
            const std::size_t br = name.find('[');
            if (br != std::string::npos)
                name = name.substr(0, br);
            if (name == "gmsh:physical")
                name += "_data";  // never shadow the property-id array
            const auto& id_map = elemental ? element_ids : condition_ids;
            StagedData sd;
            sd.mName = name;
            const int nc = mdpa_parse_data_block(
                cur, end_token, /*nodal=*/false,
                [&](std::int64_t id, std::size_t& block, std::size_t& row) {
                    auto it = id_map.find(id);
                    if (it == id_map.end())
                        return false;
                    block = it->second.first;
                    row = it->second.second;
                    return true;
                },
                sd.mRows, sd.mHasFixed);
            if (nc < 0) {
                log::warn("mdpa: {} block '{}' held no usable line", end_token.substr(4), name);
                continue;
            }
            sd.mComponents = nc;
            cell_data.push_back(std::move(sd));
        } else if (mdpa_starts_with(line, "Begin SubModelPartData")) {
            mdpa_expect_empty_block(cur, "End SubModelPartData",
                                    "a non-empty SubModelPartData block");
        } else if (mdpa_starts_with(line, "Begin SubModelPartTables")) {
            mdpa_expect_empty_block(cur, "End SubModelPartTables",
                                    "a non-empty SubModelPartTables block");
        } else if (mdpa_starts_with(line, "Begin SubModelPartGeometries")) {
            mdpa_expect_empty_block(cur, "End SubModelPartGeometries",
                                    "a non-empty SubModelPartGeometries block");
        } else if (mdpa_starts_with(line, "Begin SubModelPartConstraints")) {
            mdpa_expect_empty_block(cur, "End SubModelPartConstraints",
                                    "a non-empty SubModelPartConstraints block");
        } else if (mdpa_starts_with(line, "Begin SubModelPartNodes")) {
            if (smp_stack.empty())
                throw ReadError("MDPA: SubModelPartNodes outside a SubModelPart");
            std::vector<std::int64_t> ids;
            read_id_list("End SubModelPartNodes", ids);
            auto& smp = smps[smp_name()];
            for (std::int64_t id : ids)
                if (id >= 1 && static_cast<std::size_t>(id) <= num_points)
                    smp.mNodes.push_back(id - 1);
        } else if (mdpa_starts_with(line, "Begin SubModelPartElements") ||
                   mdpa_starts_with(line, "Begin SubModelPartConditions")) {
            const bool elemental = mdpa_starts_with(line, "Begin SubModelPartElements");
            if (smp_stack.empty())
                throw ReadError("MDPA: SubModelPart entity list outside a SubModelPart");
            std::vector<std::int64_t> ids;
            read_id_list(elemental ? "End SubModelPartElements" : "End SubModelPartConditions",
                         ids);
            const auto& id_map = elemental ? element_ids : condition_ids;
            auto& smp = smps[smp_name()];
            for (std::int64_t id : ids) {
                auto it = id_map.find(id);
                if (it == id_map.end()) {
                    log::warn("mdpa: SubModelPart references unknown entity id {}", id);
                    continue;
                }
                smp.mCells.push_back(it->second);
            }
        } else if (mdpa_starts_with(line, "Begin SubModelPart")) {
            const std::vector<std::string> head = mdpa_tokens(line);
            if (head.size() < 3)
                throw ReadError("MDPA: malformed SubModelPart header: " + line);
            smp_stack.push_back(head[2]);
            smps[smp_name()];  // an entity-less SubModelPart is still a group
        } else if (line == "End SubModelPart") {
            if (smp_stack.empty())
                throw ReadError("MDPA: 'End SubModelPart' without a matching 'Begin'");
            smp_stack.pop_back();
        } else if (mdpa_starts_with(line, "Begin Table")) {
            throw ReadError("MDPA: Table blocks are not supported by the C++ reader");
        } else if (mdpa_starts_with(line, "Begin Geometries")) {
            throw ReadError("MDPA: Geometries blocks are not supported by the C++ reader");
        } else if (mdpa_starts_with(line, "Begin Mesh")) {
            throw ReadError("MDPA: Mesh blocks are not supported by the C++ reader");
        } else if (mdpa_starts_with(line, "Begin ")) {
            throw ReadError("MDPA: unsupported block '" + line + "'");
        } else {
            throw ReadError("MDPA: unexpected line outside a block: '" + line + "'");
        }
    }
    if (!smp_stack.empty())
        throw ReadError("MDPA: EOF before 'End SubModelPart'");

    // ---- materialize ------------------------------------------------------
    Mesh mesh;
    {
        NDArray pts(DType::Float64, {num_points, 3});
        double* pp = pts.As<double>();
        for (std::size_t i = 0; i < coords.size(); ++i)
            pp[i] = coords[i];
        mesh.AssignPoints(std::move(pts));
    }

    std::vector<NDArray> props;
    std::vector<std::size_t> block_base(blocks.size(), 0);
    std::size_t running = 0;
    for (std::size_t b = 0; b < blocks.size(); ++b) {
        MdpaBlock& blk = blocks[b];
        block_base[b] = running;
        running += blk.mCount;
        for (std::int64_t node : blk.mConn)
            if (node < 0 || static_cast<std::size_t>(node) >= num_points)
                throw ReadError("MDPA: connectivity refers to node " + std::to_string(node + 1) +
                                " but the file has " + std::to_string(num_points) + " nodes");
        NDArray conn(DType::Int64, {blk.mCount, blk.mNodes});
        std::int64_t* cp = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < blk.mConn.size(); ++i)
            cp[i] = blk.mConn[i];
        mesh.AddCellBlock(blk.mType, std::move(conn));
        NDArray tag(DType::Int64, {blk.mCount});
        std::int64_t* tp = tag.As<std::int64_t>();
        for (std::size_t i = 0; i < blk.mProps.size(); ++i)
            tp[i] = blk.mProps[i];
        props.push_back(std::move(tag));
    }
    if (!blocks.empty())
        mesh.AddCellData("gmsh:physical", std::move(props));

    for (auto& fd : field_data)
        mesh.AddFieldData(fd.first, std::move(fd.second));

    for (const auto& sd : point_data) {
        mesh.AddPointData(sd.mName, mdpa_data_array(sd.mRows, 0, num_points, sd.mComponents));
        if (sd.mHasFixed) {
            NDArray fx(DType::Int64, {num_points});
            std::int64_t* fp = fx.As<std::int64_t>();
            for (std::size_t i = 0; i < num_points; ++i)
                fp[i] = -1;
            for (const auto& r : sd.mRows)
                if (r.mFixed >= 0 && r.mRow < num_points)
                    fp[r.mRow] = r.mFixed;
            mesh.AddPointData(sd.mName + "_fixed_status", std::move(fx));
        }
    }
    for (const auto& sd : cell_data) {
        std::vector<NDArray> arrays;
        arrays.reserve(blocks.size());
        for (std::size_t b = 0; b < blocks.size(); ++b)
            arrays.push_back(mdpa_data_array(sd.mRows, b, blocks[b].mCount, sd.mComponents));
        mesh.AddCellData(sd.mName, std::move(arrays));
    }

    for (const auto& smp : smps) {
        if (smp.second.mNodes.empty() && smp.second.mCells.empty()) {
            // An entity-less SubModelPart is still a named group: carry it as
            // an empty Point region rather than losing the name.
            mesh.AddRegion(Region(smp.first, RegionKind::Point, NDArray(DType::Int64, {0})));
            continue;
        }
        if (!smp.second.mNodes.empty()) {
            NDArray e(DType::Int64, {smp.second.mNodes.size()});
            std::int64_t* p = e.As<std::int64_t>();
            for (std::size_t i = 0; i < smp.second.mNodes.size(); ++i)
                p[i] = smp.second.mNodes[i];
            mesh.AddRegion(Region(smp.first, RegionKind::Point, std::move(e)));
        }
        if (!smp.second.mCells.empty()) {
            NDArray e(DType::Int64, {smp.second.mCells.size()});
            std::int64_t* p = e.As<std::int64_t>();
            for (std::size_t i = 0; i < smp.second.mCells.size(); ++i)
                p[i] = static_cast<std::int64_t>(block_base[smp.second.mCells[i].first] +
                                                 smp.second.mCells[i].second);
            mesh.AddRegion(Region(smp.first, RegionKind::Cell, std::move(e)));
        }
    }
    return mesh;
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

namespace {

bool mdpa_is_int_dtype(DType dt) {
    return dt == DType::Int8 || dt == DType::Int16 || dt == DType::Int32 || dt == DType::Int64 ||
           dt == DType::UInt8 || dt == DType::UInt16 || dt == DType::UInt32 || dt == DType::UInt64;
}

/// One value of @p rArray, formatted the way MDPA spells numbers.
std::string mdpa_format_value(const NDArray& rArray, std::size_t index) {
    char buf[64];
    if (mdpa_is_int_dtype(rArray.Dtype())) {
        std::snprintf(buf, sizeof(buf), "%lld",
                      static_cast<long long>(detail::read_int(rArray, index)));
    } else {
        std::snprintf(buf, sizeof(buf), "%.16g", detail::read_double(rArray, index));
    }
    return buf;
}

/// Number of trailing components of a data array (1 for a 1-D array).
std::size_t mdpa_components(const NDArray& rArray, std::size_t rows) {
    if (rows == 0)
        return 1;
    return rArray.Size() / rows;
}

bool mdpa_skip_point_data(const std::string& rName) {
    const std::string suffix = "_fixed_status";
    if (rName.size() >= suffix.size() &&
        rName.compare(rName.size() - suffix.size(), suffix.size(), suffix) == 0)
        return true;
    return mdpa_starts_with(rName, "gmsh:");
}

bool mdpa_skip_cell_data(const std::string& rName) {
    const std::string suffix = "_tag";
    if (rName.size() >= suffix.size() &&
        rName.compare(rName.size() - suffix.size(), suffix.size(), suffix) == 0)
        return true;
    return mdpa_starts_with(rName, "gmsh:");
}

}  // namespace

void write_mdpa(const std::string& rPath, const Mesh& rMesh) {
    std::ofstream os(rPath);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    // ---- per-block decisions (entity kind, Kratos name, written ids) -------
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<bool> is_condition(nblocks, false);
    std::vector<std::string> entity_name(nblocks);
    std::vector<std::vector<std::int64_t>> written_ids(nblocks);
    std::vector<std::size_t> block_base(nblocks, 0);
    {
        std::int64_t next_element = 1, next_condition = 1;
        std::size_t running = 0;
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            if (cb.IsRagged())
                throw WriteError("MDPA: ragged/polyhedron cell blocks are not supported (block " +
                                 std::to_string(b) + ", type '" + std::string(cb.Type()) + "')");
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            if (type == CellType::Custom)
                throw WriteError("MDPA: no Kratos entity name for cell type '" +
                                 std::string(cb.Type()) + "'");
            // The Python reference's rule for a mesh with no physical tags: a
            // block whose default Kratos *element* name is a 2-D one is written
            // as a Condition, everything else as an Element.
            is_condition[b] = kratos_element_name(type).find("2D") != std::string::npos;
            entity_name[b] =
                is_condition[b] ? kratos_condition_name(type) : kratos_element_name(type);
            block_base[b] = running;
            running += cb.NumCells();
            written_ids[b].resize(cb.NumCells());
            for (std::size_t r = 0; r < cb.NumCells(); ++r)
                written_ids[b][r] = is_condition[b] ? next_condition++ : next_element++;
        }
    }

    // ---- ModelPartData ----------------------------------------------------
    os << "Begin ModelPartData\n";
    for (const auto& name : rMesh.FieldDataNames()) {
        const NDArray& a = rMesh.FieldData(name);
        if (a.Size() != 1) {
            log::warn("mdpa: field_data '{}' has {} values; only scalars are written", name,
                      a.Size());
            continue;
        }
        os << "    " << name << " " << mdpa_format_value(a, 0) << "\n";
    }
    os << "End ModelPartData\n\n";
    os << "Begin Properties 0\nEnd Properties\n\n";

    // ---- Nodes ------------------------------------------------------------
    os << "Begin Nodes\n";
    {
        const NDArray& points = rMesh.Points();
        const std::size_t dim = rMesh.PointDim();
        const std::size_t np = rMesh.NumPoints();
        char buf[64];
        for (std::size_t i = 0; i < np; ++i) {
            os << " " << (i + 1);
            for (std::size_t c = 0; c < 3; ++c) {
                const double v = c < dim ? detail::read_double(points, i * dim + c) : 0.0;
                std::snprintf(buf, sizeof(buf), "%.16e", v);
                os << " " << buf;
            }
            os << "\n";
        }
    }
    os << "End Nodes\n\n";

    // ---- Elements / Conditions -------------------------------------------
    const bool has_props =
        rMesh.HasCellData("gmsh:physical") && rMesh.CellDataNumBlocks("gmsh:physical") == nblocks;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const CellType type = cell_type_from_name(std::string(cb.Type()));
        const std::vector<int>& order = mdpa_kratos_node_order(type);
        const std::string kind = is_condition[b] ? "Conditions" : "Elements";
        os << "Begin " << kind << " " << entity_name[b] << "\n";
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::int64_t prop = 0;
            if (has_props)
                prop = detail::read_int(rMesh.CellData("gmsh:physical", b), r);
            os << "  " << written_ids[b][r] << " " << prop;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t slot = order.empty() ? j : static_cast<std::size_t>(order[j]);
                os << " " << (detail::read_int(conn, r * k + slot) + 1);
            }
            os << "\n";
        }
        os << "End " << kind << "\n\n";
    }

    // ---- NodalData --------------------------------------------------------
    const std::size_t np = rMesh.NumPoints();
    for (const auto& name : rMesh.PointDataNames()) {
        if (mdpa_skip_point_data(name))
            continue;
        const NDArray& a = rMesh.PointData(name);
        const std::size_t nc = mdpa_components(a, np);
        const std::string fixed_name = name + "_fixed_status";
        const bool has_fixed = rMesh.HasPointData(fixed_name);
        os << "Begin NodalData " << name << "\n";
        for (std::size_t i = 0; i < np; ++i) {
            bool all_nan = true;
            for (std::size_t j = 0; j < nc; ++j)
                if (!std::isnan(detail::read_double(a, i * nc + j)))
                    all_nan = false;
            if (all_nan)
                continue;
            os << "  " << (i + 1);
            if (has_fixed) {
                const std::int64_t f = detail::read_int(rMesh.PointData(fixed_name), i);
                if (f >= 0)
                    os << " " << f;
            }
            for (std::size_t j = 0; j < nc; ++j)
                os << " " << mdpa_format_value(a, i * nc + j);
            os << "\n";
        }
        os << "End NodalData\n\n";
    }

    // ---- ElementalData / ConditionalData ----------------------------------
    for (const auto& name : rMesh.CellDataNames()) {
        if (mdpa_skip_cell_data(name))
            continue;
        if (rMesh.CellDataNumBlocks(name) != nblocks)
            continue;
        for (int pass = 0; pass < 2; ++pass) {
            const bool conditions = pass == 1;
            std::ostringstream body;
            for (std::size_t b = 0; b < nblocks; ++b) {
                if (is_condition[b] != conditions)
                    continue;
                const auto cb = rMesh.Cells(b);
                const NDArray& a = rMesh.CellData(name, b);
                const std::size_t nc = mdpa_components(a, cb.NumCells());
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    bool all_nan = true;
                    for (std::size_t j = 0; j < nc; ++j)
                        if (!std::isnan(detail::read_double(a, r * nc + j)))
                            all_nan = false;
                    if (all_nan)
                        continue;
                    body << "  " << written_ids[b][r];
                    for (std::size_t j = 0; j < nc; ++j)
                        body << " " << mdpa_format_value(a, r * nc + j);
                    body << "\n";
                }
            }
            if (body.str().empty())
                continue;
            const std::string kind = conditions ? "ConditionalData" : "ElementalData";
            os << "Begin " << kind << " " << name << "\n" << body.str() << "End " << kind << "\n\n";
        }
    }

    // ---- SubModelParts from named regions ---------------------------------
    for (const auto& name : rMesh.RegionNames()) {
        std::vector<std::int64_t> nodes;
        std::vector<std::int64_t> elements, conditions;
        bool any = false;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = rMesh.Region(i);
            if (r.mName != name)
                continue;
            if (r.mKind == RegionKind::Side) {
                log::warn("mdpa: dropping side region '{}' (MDPA has no facet sets)", name);
                continue;
            }
            any = true;
            const std::int64_t* e = r.Entries();
            for (std::size_t j = 0; j < r.NumEntries(); ++j) {
                if (r.mKind == RegionKind::Point) {
                    nodes.push_back(e[j] + 1);
                    continue;
                }
                // Cell region: global block-major index -> (block, row) -> id.
                const std::size_t g = static_cast<std::size_t>(e[j]);
                for (std::size_t b = 0; b < nblocks; ++b) {
                    const std::size_t count = rMesh.Cells(b).NumCells();
                    if (g < block_base[b] || g >= block_base[b] + count)
                        continue;
                    const std::int64_t id = written_ids[b][g - block_base[b]];
                    (is_condition[b] ? conditions : elements).push_back(id);
                    break;
                }
            }
        }
        if (!any)
            continue;
        os << "Begin SubModelPart " << name << "\n";
        auto emit = [&](const char* pTag, const std::vector<std::int64_t>& rIds) {
            if (rIds.empty())
                return;
            os << "    Begin SubModelPart" << pTag << "\n";
            for (std::int64_t id : rIds)
                os << "        " << id << "\n";
            os << "    End SubModelPart" << pTag << "\n";
        };
        emit("Nodes", nodes);
        emit("Elements", elements);
        emit("Conditions", conditions);
        os << "End SubModelPart\n\n";
    }
}

}  // namespace meshioplusplus
