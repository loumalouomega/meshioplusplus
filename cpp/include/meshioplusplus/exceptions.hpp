#pragma once
//
// meshio I/O exceptions. The binding layer translates these to the existing
// Python meshio.ReadError / meshio.WriteError classes so behaviour is
// identical from Python.

#include <stdexcept>
#include <string>

namespace meshioplusplus {

struct ReadError : std::runtime_error {
    ReadError() : std::runtime_error("") {}
    explicit ReadError(const std::string& msg) : std::runtime_error(msg) {}
};

struct WriteError : std::runtime_error {
    WriteError() : std::runtime_error("") {}
    explicit WriteError(const std::string& msg) : std::runtime_error(msg) {}
};

}  // namespace meshioplusplus
