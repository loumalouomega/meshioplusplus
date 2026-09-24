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
#include <cstdint>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/binary_stream.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The marker at `Offset`, or -1 when it does not fit.
std::int64_t fortran_marker(const char* pData, std::size_t Size, std::size_t Offset,
                            const FortranRecordLayout& rLayout) {
    const std::size_t w = static_cast<std::size_t>(rLayout.mMarkerBytes);
    if (Offset > Size || Size - Offset < w)
        return -1;
    ByteCursor c(pData + Offset, w, rLayout.mBigEndian, "Fortran records");
    const std::int64_t v = w == 4 ? static_cast<std::int64_t>(c.I32()) : c.I64();
    return v;
}

}  // namespace

std::optional<FortranRecordLayout> sniff_fortran_records(const char* pData, std::size_t Size) {
    for (int bytes : {4, 8}) {
        for (bool big : {false, true}) {
            const FortranRecordLayout layout{bytes, big};
            const std::int64_t n = fortran_marker(pData, Size, 0, layout);
            if (n <= 0)
                continue;
            const std::size_t w = static_cast<std::size_t>(bytes);
            const std::size_t len = static_cast<std::size_t>(n);
            if (len > Size - w)
                continue;
            if (fortran_marker(pData, Size, w + len, layout) == n)
                return layout;
        }
    }
    return std::nullopt;
}

std::vector<FortranRecord> fortran_records(const char* pData, std::size_t Size,
                                           const FortranRecordLayout& rLayout,
                                           const std::string& rWhat) {
    std::vector<FortranRecord> out;
    const std::size_t w = static_cast<std::size_t>(rLayout.mMarkerBytes);
    std::size_t pos = 0;
    while (pos < Size) {
        const std::int64_t n = fortran_marker(pData, Size, pos, rLayout);
        if (n < 0 || static_cast<std::uint64_t>(n) > Size - pos - w)
            throw ReadError(rWhat + ": Fortran record at offset " + std::to_string(pos) +
                            " runs past the end of the file");
        const std::size_t len = static_cast<std::size_t>(n);
        if (fortran_marker(pData, Size, pos + w + len, rLayout) != n)
            throw ReadError(rWhat + ": Fortran record at offset " + std::to_string(pos) +
                            " has mismatched length markers");
        out.push_back({pos + w, len});
        pos += 2 * w + len;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
