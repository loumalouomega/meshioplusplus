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
#pragma once

/**
 * @file detail/read_guard.hpp
 * @brief Makes every reader entry point fail with ReadError, whatever the parser threw.
 *
 * The format readers parse with `std::stoll`, `std::vector::at` and friends,
 * which report a malformed token as `std::invalid_argument`/`std::out_of_range`
 * rather than as the `ReadError` the fallback contract is built on (pybind11
 * would surface them as `ValueError`/`IndexError`, and `core_declined` would
 * log them as a broken fast path). Wrapping each entry point -- the registry's
 * readers and the Python `*_read` bindings, not the global exception
 * translator, which must keep mapping an operation's `std::invalid_argument`
 * to `ValueError` -- makes "the file is not something this reader reads" one
 * exception on every surface. `ReadError`, `WriteError`, `Unsupported` and
 * `std::bad_alloc` pass through unchanged, as does anything that is not a
 * `std::exception` subclass the parser could have thrown.
 */

// System includes
#include <exception>
#include <new>
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus::detail {

/**
 * @brief Rethrows the in-flight exception, as ReadError when a parser threw it.
 *
 * Call only from inside a `catch (...)` block.
 * @param pFormat The format name, used in the message ("malformed <fmt> input").
 */
[[noreturn]] inline void rethrow_as_read_error(const char* pFormat) {
    try {
        throw;
    } catch (const ReadError&) {
        throw;
    } catch (const WriteError&) {
        throw;
    } catch (const Unsupported&) {
        throw;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::logic_error& e) {
        throw ReadError(std::string("meshio++: malformed ") + pFormat + " input (" + e.what() +
                        ")");
    } catch (const std::runtime_error& e) {
        throw ReadError(std::string("meshio++: ") + pFormat + ": " + e.what());
    }
}

/**
 * @brief Calls @p rFn, translating a parser's exceptions with rethrow_as_read_error().
 * @param pFormat The format name for the message.
 * @param rFn The read to run.
 * @return Whatever @p rFn returns.
 */
template <class TFn>
decltype(auto) guarded_read(const char* pFormat, TFn&& rFn) {
    try {
        return rFn();
    } catch (...) {
        rethrow_as_read_error(pFormat);
    }
}

}  // namespace meshioplusplus::detail
