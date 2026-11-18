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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#include "integer_precision.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>

#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// Which arrays a format carries at all (an array it drops is not worth a warning).
struct IpKinds {
    bool mPoint;
    bool mCell;
};

// Keep in step with `_FLOAT_STORED_INTEGER_FORMATS` in
// src/python/meshioplusplus/_helpers.py; test_integer_precision.py compares them.
const std::map<std::string, IpKinds>& ip_formats() {
    static const std::map<std::string, IpKinds> m = {
        {"cae", {false, true}},          {"ensight", {true, true}},
        {"frd", {true, false}},          {"gid", {true, true}},
        {"gltf", {true, false}},         {"marc_t19", {true, true}},
        {"radioss_anim", {true, true}},  {"xplt", {true, true}},
    };
    return m;
}

bool ip_inexact(const NDArray& rA) {
    if (is_float_dtype(rA.Dtype()))
        return false;
    const std::size_t n = rA.Size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::int64_t v = read_int(rA, i);
        if (v > kExactFloatInteger || v < -kExactFloatInteger)
            return true;
    }
    return false;
}

void ip_warn(const std::string& rFormat, const std::string& rName) {
    log::warn("{}: the integer array '{}' is stored as floating point, so values above {} are rounded",
              rFormat, rName, kExactFloatInteger);
}

}  // namespace

void warn_integer_precision(const std::string& rFormat, const Mesh& rMesh) {
    const auto it = ip_formats().find(rFormat);
    if (it == ip_formats().end())
        return;
    if (it->second.mPoint) {
        for (const std::string& name : rMesh.PointDataNames())
            if (ip_inexact(rMesh.PointData(name)))
                ip_warn(rFormat, name);
    }
    if (it->second.mCell) {
        for (const std::string& name : rMesh.CellDataNames()) {
            for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                if (ip_inexact(rMesh.CellData(name, b))) {
                    ip_warn(rFormat, name);
                    break;
                }
            }
        }
    }
}

}  // namespace detail
}  // namespace meshioplusplus
