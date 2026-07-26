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
 * @file operations/data_calc.hpp
 * @brief Derives a new data array from an elementwise expression over the
 * existing arrays at one location.
 *
 * The evaluator is a small hand-written tokenizer plus recursive-descent
 * parser — there is **no** external parser library and **no** evaluation of
 * arbitrary code. Only the tokens below are accepted; anything else is a
 * diagnosed error.
 *
 * ### Grammar
 * @code
 * expr    := term (('+'|'-') term)*
 * term    := unary (('*'|'/') unary)*
 * unary   := ('-'|'+') unary | primary
 * primary := number | ident | ident '(' expr (',' expr)* ')' | '(' expr ')'
 * @endcode
 *
 * Functions: `abs(x)`, `sqrt(x)`, `min(a,b)`, `max(a,b)`, and `norm(v)` for the
 * Euclidean magnitude of a vector. Exponentiation is deliberately absent.
 *
 * Identifiers may contain `:` and `.` after their first character, so a mesh's
 * own conventional names (`gmsh:physical`, `quality:scaled_jacobian`) can be
 * referenced directly. A name containing spaces or operator characters can be
 * written between backticks: `` `my array` ``.
 *
 * ### Values
 * Every operand is either a scalar (1 component) or a vector/tensor of up to
 * `DATA_CALC_MAX_COMPONENTS` components. Binary operators require matching
 * component counts, except that a scalar broadcasts against a wider operand.
 * `norm` always produces a scalar. Arithmetic is performed in `double` and the
 * result is stored in `DataCalcOptions::dtype`.
 *
 * Division by zero yields an IEEE infinity or NaN and is not an error, per the
 * policy in `operations/data_common.hpp`.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format — it is not
 * in the format registry.
 */

// System includes
#include <cstddef>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {

/// Widest value the evaluator handles — covers scalars, 3-vectors and 3x3
/// tensors. Bounding this is what lets the element loop use stack buffers and
/// allocate nothing.
inline constexpr std::size_t DATA_CALC_MAX_COMPONENTS = 16;

/// Deepest expression nesting accepted, guarding the recursive-descent parser
/// against a stack overflow from hostile input (expressions reach this code
/// from the C ABI and the CLIs).
inline constexpr std::size_t DATA_CALC_MAX_DEPTH = 64;

/// Options for `data_calc`.
struct DataCalcOptions {
    /// Where operands are looked up and the result is stored.
    DataLocation location = DataLocation::Point;
    /// Name of the new array. Required.
    std::string output;
    /// Whether an existing array of that name may be replaced.
    bool overwrite = false;
    /// Dtype of the produced array; the arithmetic itself is always `double`.
    DType dtype = DType::Float64;
};

/**
 * @brief Evaluates @p rExpression elementwise and stores the result as a new
 * array at `rOpts.location`.
 *
 * For `DataLocation::Cell` the expression is evaluated once per cell block, so
 * the produced `cell_data` always has exactly one array per block.
 * @param rMesh the source mesh (unmodified).
 * @param rExpression the expression text (see the file comment for the grammar).
 * @param rOpts the output name, location and dtype.
 * @return a new mesh carrying the derived array alongside everything the input
 *         already had.
 * @throws std::invalid_argument on any lexical, syntactic, name-resolution,
 *         arity, component-width or row-count error. Every message is prefixed
 *         `meshio++: data_calc: ` and, where meaningful, carries the 0-based
 *         character position within the expression.
 */
MESHIOPLUSPLUS_API Mesh data_calc(const Mesh& rMesh, const std::string& rExpression, const DataCalcOptions& rOpts);

}  // namespace meshioplusplus
