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
// Shared vocabulary for the data operations: the DataLocation / NanPolicy
// string parsers used by both CLIs and the WASM binding, the sorted name
// listing, and the "unknown key (available: ...)" diagnostic every data
// operation reports. See operations/data_common.hpp for the contract.

// System includes
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

DataLocation data_location_from_name(const std::string& rName) {
    if (rName == "point" || rName == "point_data" || rName == "node")
        return DataLocation::Point;
    if (rName == "cell" || rName == "cell_data" || rName == "element")
        return DataLocation::Cell;
    if (rName == "field" || rName == "field_data")
        return DataLocation::Field;
    throw std::invalid_argument("meshio++: unknown data location '" + rName +
                                "' (expected 'point', 'cell' or 'field')");
}

const char* data_location_name(DataLocation Location) {
    switch (Location) {
        case DataLocation::Point:
            return "point_data";
        case DataLocation::Cell:
            return "cell_data";
        case DataLocation::Field:
            return "field_data";
    }
    return "point_data";  // unreachable
}

NanPolicy nan_policy_from_name(const std::string& rName) {
    if (rName == "ignore")
        return NanPolicy::Ignore;
    if (rName == "replace")
        return NanPolicy::Replace;
    if (rName == "fail")
        return NanPolicy::Fail;
    throw std::invalid_argument("meshio++: unknown NaN policy '" + rName +
                                "' (expected 'ignore', 'replace' or 'fail')");
}

std::vector<std::string> data_names(const Mesh& rMesh, DataLocation Location) {
    switch (Location) {
        case DataLocation::Point:
            return rMesh.PointDataNames();
        case DataLocation::Cell:
            return rMesh.CellDataNames();
        case DataLocation::Field:
            return rMesh.FieldDataNames();
    }
    return {};  // unreachable
}

bool data_has(const Mesh& rMesh, DataLocation Location, const std::string& rName) {
    switch (Location) {
        case DataLocation::Point:
            return rMesh.HasPointData(rName);
        case DataLocation::Cell:
            return rMesh.HasCellData(rName);
        case DataLocation::Field:
            return rMesh.HasFieldData(rName);
    }
    return false;  // unreachable
}

std::string data_unknown_key_message(const Mesh& rMesh, DataLocation Location,
                                     const std::string& rName) {
    const char* loc = data_location_name(Location);
    const std::vector<std::string> available = data_names(rMesh, Location);
    std::string msg = "meshio++: no ";
    msg += loc;
    msg += " array named '" + rName + "' (";
    if (available.empty()) {
        msg += "the mesh has no ";
        msg += loc;
    } else {
        msg += "available: ";
        for (std::size_t i = 0; i < available.size(); ++i) {
            if (i != 0)
                msg += ", ";
            msg += available[i];
        }
    }
    msg += ")";
    return msg;
}

std::size_t data_num_components(const NDArray& rArray) {
    const std::vector<std::size_t>& shape = rArray.Shape();
    if (shape.size() < 2)
        return 1;
    std::size_t n = 1;
    for (std::size_t d = 1; d < shape.size(); ++d)
        n *= shape[d];
    return n == 0 ? 1 : n;
}

}  // namespace meshioplusplus
