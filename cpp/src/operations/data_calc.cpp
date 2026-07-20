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
// Elementwise expression evaluator over a mesh's data arrays. Hand-written
// tokenizer + recursive-descent parser; no external parser library, no
// evaluation of arbitrary code -- only the fixed token and function set in
// operations/data_calc.hpp is accepted.
//
// The parse produces a FLAT node pool (std::vector<CalcNode>) rather than a
// pointer tree: it is trivially shareable read-only across parallel_for
// threads, needs no per-node allocation, and keeps evaluation cache-friendly.
//
// Evaluation runs in three passes:
//   A. infer each node's component width and resolve array names (once);
//   B. widen every referenced array to double (once per array, not per node);
//   C. the element loop, which is branch-free over dtypes and allocates
//      nothing because widths are bounded by DATA_CALC_MAX_COMPONENTS.
//
// Numeric literals go through std::strtod rather than std::from_chars: the
// floating-point overload of the latter is missing from some libc++ versions
// (notably a hazard for the Emscripten/WASM build), and strtod matches the
// existing precedent in the format readers (dex.cpp, wkt.cpp, nastran.cpp).
//
// Geometry is never modified. See operations/data_calc.hpp for the contract.

// System includes
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

[[noreturn]] void calc_fail(const std::string& rMessage) {
    throw std::invalid_argument("meshio++: data_calc: " + rMessage);
}

[[noreturn]] void calc_fail_at(const std::string& rMessage, std::size_t pos) {
    calc_fail(rMessage + " at position " + std::to_string(pos));
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

enum class CalcTokenType {
    Number,
    Ident,
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
    Comma,
    End,
};

struct CalcToken {
    CalcTokenType mType = CalcTokenType::End;
    double mValue = 0.0;   ///< Set for Number.
    std::string mText;     ///< Set for Ident.
    std::size_t mPos = 0;  ///< 0-based offset of the token's first character.
};

bool calc_ident_start(char c) {
    return std::isalpha(static_cast<unsigned char>(c)) != 0 || c == '_';
}

bool calc_ident_body(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_' || c == ':' || c == '.';
}

std::vector<CalcToken> calc_tokenize(const std::string& rText) {
    std::vector<CalcToken> out;
    std::size_t i = 0;
    const std::size_t n = rText.size();
    while (i < n) {
        const char c = rText[i];
        if (std::isspace(static_cast<unsigned char>(c)) != 0) {
            ++i;
            continue;
        }
        CalcToken t;
        t.mPos = i;
        if (std::isdigit(static_cast<unsigned char>(c)) != 0 ||
            (c == '.' && i + 1 < n &&
             std::isdigit(static_cast<unsigned char>(rText[i + 1])) != 0)) {
            const char* begin = rText.c_str() + i;
            char* end = nullptr;
            const double v = std::strtod(begin, &end);
            if (end == begin)
                calc_fail_at("malformed number", i);
            t.mType = CalcTokenType::Number;
            t.mValue = v;
            i += static_cast<std::size_t>(end - begin);
            out.push_back(std::move(t));
            continue;
        }
        if (c == '`') {
            // Backtick-quoted identifier, for names with spaces or operator
            // characters that the bare identifier rule cannot express.
            const std::size_t close = rText.find('`', i + 1);
            if (close == std::string::npos)
                calc_fail_at("unterminated `-quoted name", i);
            t.mType = CalcTokenType::Ident;
            t.mText = rText.substr(i + 1, close - i - 1);
            if (t.mText.empty())
                calc_fail_at("empty `-quoted name", i);
            i = close + 1;
            out.push_back(std::move(t));
            continue;
        }
        if (calc_ident_start(c)) {
            std::size_t j = i + 1;
            while (j < n && calc_ident_body(rText[j]))
                ++j;
            // A trailing ':' or '.' is punctuation, not part of the name.
            while (j > i + 1 && (rText[j - 1] == ':' || rText[j - 1] == '.'))
                --j;
            t.mType = CalcTokenType::Ident;
            t.mText = rText.substr(i, j - i);
            i = j;
            out.push_back(std::move(t));
            continue;
        }
        switch (c) {
            case '+':
                t.mType = CalcTokenType::Plus;
                break;
            case '-':
                t.mType = CalcTokenType::Minus;
                break;
            case '*':
                t.mType = CalcTokenType::Star;
                break;
            case '/':
                t.mType = CalcTokenType::Slash;
                break;
            case '(':
                t.mType = CalcTokenType::LParen;
                break;
            case ')':
                t.mType = CalcTokenType::RParen;
                break;
            case ',':
                t.mType = CalcTokenType::Comma;
                break;
            default:
                calc_fail_at(std::string("unexpected character '") + c + "'", i);
        }
        ++i;
        out.push_back(std::move(t));
    }
    CalcToken end;
    end.mType = CalcTokenType::End;
    end.mPos = n;
    out.push_back(std::move(end));
    return out;
}

// ---------------------------------------------------------------------------
// Parser -> flat node pool
// ---------------------------------------------------------------------------

enum class CalcNodeType {
    Const,
    Array,
    Neg,
    Add,
    Sub,
    Mul,
    Div,
    Abs,
    Sqrt,
    Min,
    Max,
    Norm,
};

struct CalcNode {
    CalcNodeType mType = CalcNodeType::Const;
    double mConst = 0.0;
    std::size_t mOperand = 0;  ///< Index into the resolved-operand table (Array).
    std::size_t mLhs = 0;      ///< Index into the node pool.
    std::size_t mRhs = 0;      ///< Index into the node pool.
    std::size_t mWidth = 1;    ///< Component count, filled by the inference pass.
    std::size_t mPos = 0;      ///< Source position, for diagnostics.
    std::string mName;         ///< Array name (Array nodes only).
};

/// The known function names, and their arity.
struct CalcFuncDef {
    const char* mName;
    CalcNodeType mType;
    int mArity;
};

const CalcFuncDef CALC_FUNCS[] = {
    {"abs", CalcNodeType::Abs, 1}, {"sqrt", CalcNodeType::Sqrt, 1}, {"min", CalcNodeType::Min, 2},
    {"max", CalcNodeType::Max, 2}, {"norm", CalcNodeType::Norm, 1},
};

std::string calc_known_funcs() {
    std::string s;
    for (const CalcFuncDef& f : CALC_FUNCS) {
        if (!s.empty())
            s += ", ";
        s += f.mName;
    }
    return s;
}

class CalcParser {
public:
    CalcParser(const std::vector<CalcToken>& rTokens, std::vector<CalcNode>& rPool)
        : mrTokens(rTokens), mrPool(rPool) {}

    /// Parses the whole token stream; returns the root node index.
    std::size_t ParseAll() {
        const std::size_t root = ParseExpr(0);
        if (Peek().mType != CalcTokenType::End)
            calc_fail_at("trailing input after the expression", Peek().mPos);
        return root;
    }

private:
    const CalcToken& Peek() const { return mrTokens[mPos]; }
    const CalcToken& Take() { return mrTokens[mPos++]; }

    std::size_t Emit(CalcNode node) {
        mrPool.push_back(std::move(node));
        return mrPool.size() - 1;
    }

    void CheckDepth(std::size_t depth) {
        if (depth > DATA_CALC_MAX_DEPTH)
            calc_fail("expression nests deeper than " + std::to_string(DATA_CALC_MAX_DEPTH) +
                      " levels");
    }

    std::size_t ParseExpr(std::size_t depth) {
        CheckDepth(depth);
        std::size_t lhs = ParseTerm(depth + 1);
        for (;;) {
            const CalcTokenType t = Peek().mType;
            if (t != CalcTokenType::Plus && t != CalcTokenType::Minus)
                return lhs;
            const CalcToken& op = Take();
            const std::size_t rhs = ParseTerm(depth + 1);
            CalcNode n;
            n.mType = t == CalcTokenType::Plus ? CalcNodeType::Add : CalcNodeType::Sub;
            n.mLhs = lhs;
            n.mRhs = rhs;
            n.mPos = op.mPos;
            lhs = Emit(std::move(n));
        }
    }

    std::size_t ParseTerm(std::size_t depth) {
        CheckDepth(depth);
        std::size_t lhs = ParseUnary(depth + 1);
        for (;;) {
            const CalcTokenType t = Peek().mType;
            if (t != CalcTokenType::Star && t != CalcTokenType::Slash)
                return lhs;
            const CalcToken& op = Take();
            const std::size_t rhs = ParseUnary(depth + 1);
            CalcNode n;
            n.mType = t == CalcTokenType::Star ? CalcNodeType::Mul : CalcNodeType::Div;
            n.mLhs = lhs;
            n.mRhs = rhs;
            n.mPos = op.mPos;
            lhs = Emit(std::move(n));
        }
    }

    std::size_t ParseUnary(std::size_t depth) {
        CheckDepth(depth);
        const CalcTokenType t = Peek().mType;
        if (t == CalcTokenType::Plus) {
            Take();
            return ParseUnary(depth + 1);
        }
        if (t == CalcTokenType::Minus) {
            const CalcToken& op = Take();
            const std::size_t operand = ParseUnary(depth + 1);
            CalcNode n;
            n.mType = CalcNodeType::Neg;
            n.mLhs = operand;
            n.mPos = op.mPos;
            return Emit(std::move(n));
        }
        return ParsePrimary(depth + 1);
    }

    std::size_t ParsePrimary(std::size_t depth) {
        CheckDepth(depth);
        const CalcToken& tok = Peek();
        if (tok.mType == CalcTokenType::End)
            calc_fail("unexpected end of expression");
        if (tok.mType == CalcTokenType::Number) {
            Take();
            CalcNode n;
            n.mType = CalcNodeType::Const;
            n.mConst = tok.mValue;
            n.mPos = tok.mPos;
            return Emit(std::move(n));
        }
        if (tok.mType == CalcTokenType::LParen) {
            Take();
            const std::size_t inner = ParseExpr(depth + 1);
            if (Peek().mType != CalcTokenType::RParen)
                calc_fail_at("expected ')'", Peek().mPos);
            Take();
            return inner;
        }
        if (tok.mType == CalcTokenType::Ident) {
            const CalcToken ident = Take();
            if (Peek().mType != CalcTokenType::LParen) {
                CalcNode n;
                n.mType = CalcNodeType::Array;
                n.mName = ident.mText;
                n.mPos = ident.mPos;
                return Emit(std::move(n));
            }
            // A call: the name must be one of the fixed functions.
            const CalcFuncDef* def = nullptr;
            for (const CalcFuncDef& f : CALC_FUNCS)
                if (ident.mText == f.mName)
                    def = &f;
            if (def == nullptr)
                calc_fail_at(
                    "unknown function '" + ident.mText + "' (known: " + calc_known_funcs() + ")",
                    ident.mPos);
            Take();  // '('
            std::vector<std::size_t> args;
            if (Peek().mType != CalcTokenType::RParen) {
                for (;;) {
                    args.push_back(ParseExpr(depth + 1));
                    if (Peek().mType != CalcTokenType::Comma)
                        break;
                    Take();
                }
            }
            if (Peek().mType != CalcTokenType::RParen)
                calc_fail_at("expected ')'", Peek().mPos);
            Take();
            if (static_cast<int>(args.size()) != def->mArity)
                calc_fail_at("'" + ident.mText + "' takes exactly " + std::to_string(def->mArity) +
                                 (def->mArity == 1 ? " argument (got " : " arguments (got ") +
                                 std::to_string(args.size()) + ")",
                             ident.mPos);
            CalcNode n;
            n.mType = def->mType;
            n.mLhs = args[0];
            if (args.size() > 1)
                n.mRhs = args[1];
            n.mPos = ident.mPos;
            return Emit(std::move(n));
        }
        calc_fail_at("unexpected token", tok.mPos);
    }

    const std::vector<CalcToken>& mrTokens;
    std::vector<CalcNode>& mrPool;
    std::size_t mPos = 0;
};

// ---------------------------------------------------------------------------
// Pass A: width inference + operand resolution
// ---------------------------------------------------------------------------

/// One distinct array referenced by the expression.
struct CalcOperand {
    std::string mName;
    std::size_t mWidth = 1;
    /// Values widened to double, laid out row-major as row * mWidth + k.
    /// Filled per cell block for the Cell location.
    std::vector<double> mValues;
};

const char* calc_node_op_name(CalcNodeType t) {
    switch (t) {
        case CalcNodeType::Add:
            return "+";
        case CalcNodeType::Sub:
            return "-";
        case CalcNodeType::Mul:
            return "*";
        case CalcNodeType::Div:
            return "/";
        case CalcNodeType::Min:
            return "min";
        case CalcNodeType::Max:
            return "max";
        default:
            return "?";
    }
}

/// Combines two operand widths, allowing a scalar to broadcast.
std::size_t calc_combine_width(std::size_t a, std::size_t b, CalcNodeType t, std::size_t pos) {
    if (a == b)
        return a;
    if (a == 1)
        return b;
    if (b == 1)
        return a;
    calc_fail_at("cannot combine a " + std::to_string(a) + "-component array with a " +
                     std::to_string(b) + "-component array in '" + calc_node_op_name(t) + "'",
                 pos);
}

/// Walks the pool, resolving array names and filling in every node's width.
/// @return the widths of the distinct operands, in `rOperands` order.
void calc_infer(std::vector<CalcNode>& rPool, std::size_t root, const Mesh& rMesh,
                DataLocation location, std::vector<CalcOperand>& rOperands) {
    std::unordered_map<std::string, std::size_t> seen;
    // The pool is built bottom-up by the parser, so a single forward sweep
    // visits every node after its children.
    for (std::size_t i = 0; i <= root; ++i) {
        CalcNode& n = rPool[i];
        switch (n.mType) {
            case CalcNodeType::Const:
                n.mWidth = 1;
                break;
            case CalcNodeType::Array: {
                if (!data_has(rMesh, location, n.mName))
                    calc_fail_at(
                        "unknown " + std::string(data_location_name(location)) + " array '" +
                            n.mName + "' (" +
                            [&] {
                                const std::vector<std::string> av = data_names(rMesh, location);
                                if (av.empty())
                                    return std::string("the mesh has no ") +
                                           data_location_name(location);
                                std::string s = "available: ";
                                for (std::size_t k = 0; k < av.size(); ++k) {
                                    if (k != 0)
                                        s += ", ";
                                    s += av[k];
                                }
                                return s;
                            }() +
                            ")",
                        n.mPos);
                const auto it = seen.find(n.mName);
                if (it != seen.end()) {
                    n.mOperand = it->second;
                    n.mWidth = rOperands[it->second].mWidth;
                    break;
                }
                const NDArray& arr = location == DataLocation::Point  ? rMesh.PointData(n.mName)
                                     : location == DataLocation::Cell ? rMesh.CellData(n.mName, 0)
                                                                      : rMesh.FieldData(n.mName);
                CalcOperand op;
                op.mName = n.mName;
                op.mWidth = data_num_components(arr);
                if (op.mWidth > DATA_CALC_MAX_COMPONENTS)
                    calc_fail_at("array '" + n.mName + "' has " + std::to_string(op.mWidth) +
                                     " components, more than the supported maximum of " +
                                     std::to_string(DATA_CALC_MAX_COMPONENTS),
                                 n.mPos);
                n.mOperand = rOperands.size();
                n.mWidth = op.mWidth;
                seen.emplace(n.mName, rOperands.size());
                rOperands.push_back(std::move(op));
                break;
            }
            case CalcNodeType::Neg:
            case CalcNodeType::Abs:
            case CalcNodeType::Sqrt:
                n.mWidth = rPool[n.mLhs].mWidth;
                break;
            case CalcNodeType::Norm:
                n.mWidth = 1;
                break;
            case CalcNodeType::Add:
            case CalcNodeType::Sub:
            case CalcNodeType::Mul:
            case CalcNodeType::Div:
            case CalcNodeType::Min:
            case CalcNodeType::Max:
                n.mWidth =
                    calc_combine_width(rPool[n.mLhs].mWidth, rPool[n.mRhs].mWidth, n.mType, n.mPos);
                break;
        }
    }
}

// ---------------------------------------------------------------------------
// Pass C: the element evaluator
// ---------------------------------------------------------------------------

/// Evaluates node `idx` for row `row`, writing up to DATA_CALC_MAX_COMPONENTS
/// values into `pOut`. Returns the width actually produced.
std::size_t calc_eval(const std::vector<CalcNode>& rPool, std::size_t idx,
                      const std::vector<CalcOperand>& rOperands, std::size_t row, double* pOut) {
    const CalcNode& n = rPool[idx];
    switch (n.mType) {
        case CalcNodeType::Const:
            pOut[0] = n.mConst;
            return 1;
        case CalcNodeType::Array: {
            const CalcOperand& op = rOperands[n.mOperand];
            const double* src = op.mValues.data() + row * op.mWidth;
            for (std::size_t k = 0; k < op.mWidth; ++k)
                pOut[k] = src[k];
            return op.mWidth;
        }
        case CalcNodeType::Neg: {
            const std::size_t w = calc_eval(rPool, n.mLhs, rOperands, row, pOut);
            for (std::size_t k = 0; k < w; ++k)
                pOut[k] = -pOut[k];
            return w;
        }
        case CalcNodeType::Abs: {
            const std::size_t w = calc_eval(rPool, n.mLhs, rOperands, row, pOut);
            for (std::size_t k = 0; k < w; ++k)
                pOut[k] = std::fabs(pOut[k]);
            return w;
        }
        case CalcNodeType::Sqrt: {
            const std::size_t w = calc_eval(rPool, n.mLhs, rOperands, row, pOut);
            for (std::size_t k = 0; k < w; ++k)
                pOut[k] = std::sqrt(pOut[k]);
            return w;
        }
        case CalcNodeType::Norm: {
            double buf[DATA_CALC_MAX_COMPONENTS];
            const std::size_t w = calc_eval(rPool, n.mLhs, rOperands, row, buf);
            double acc = 0.0;
            for (std::size_t k = 0; k < w; ++k)
                acc += buf[k] * buf[k];
            pOut[0] = std::sqrt(acc);
            return 1;
        }
        default:
            break;
    }
    // Binary operators, with scalar broadcast on either side.
    double lhs[DATA_CALC_MAX_COMPONENTS];
    double rhs[DATA_CALC_MAX_COMPONENTS];
    const std::size_t lw = calc_eval(rPool, n.mLhs, rOperands, row, lhs);
    const std::size_t rw = calc_eval(rPool, n.mRhs, rOperands, row, rhs);
    const std::size_t w = lw > rw ? lw : rw;
    for (std::size_t k = 0; k < w; ++k) {
        const double a = lhs[lw == 1 ? 0 : k];
        const double b = rhs[rw == 1 ? 0 : k];
        switch (n.mType) {
            case CalcNodeType::Add:
                pOut[k] = a + b;
                break;
            case CalcNodeType::Sub:
                pOut[k] = a - b;
                break;
            case CalcNodeType::Mul:
                pOut[k] = a * b;
                break;
            case CalcNodeType::Div:
                pOut[k] = a / b;
                break;
            case CalcNodeType::Min:
                pOut[k] = a < b ? a : b;
                break;
            case CalcNodeType::Max:
                pOut[k] = a > b ? a : b;
                break;
            default:
                pOut[k] = 0.0;
                break;
        }
    }
    return w;
}

/// Pass B for one group of rows: widen every referenced array to double.
/// @param rSource supplies the NDArray for an operand (differs per cell block).
template <class TSource>
void calc_load_operands(std::vector<CalcOperand>& rOperands, std::size_t rows, TSource&& rSource,
                        const char* pWhat) {
    for (CalcOperand& op : rOperands) {
        const NDArray& arr = rSource(op.mName);
        const std::size_t got = detail::rows(arr);
        if (got != rows)
            calc_fail("array '" + op.mName + "' has " + std::to_string(got) + " rows but " + pWhat +
                      " needs " + std::to_string(rows));
        if (data_num_components(arr) != op.mWidth)
            calc_fail("array '" + op.mName + "' has inconsistent component counts across blocks");
        op.mValues.resize(rows * op.mWidth);
        double* dst = op.mValues.data();
        // Bandwidth-bound widening gather, not compute.
        parallel_for_bw(rows * op.mWidth,
                        [&](std::size_t i) { dst[i] = detail::read_double(arr, i); });
    }
}

/// Runs the element loop over `rows` entries into a fresh array.
NDArray calc_run(const std::vector<CalcNode>& rPool, std::size_t root,
                 const std::vector<CalcOperand>& rOperands, std::size_t rows, std::size_t width,
                 DType dtype) {
    NDArray dst = NDArray::Uninit(
        dtype, width <= 1 ? std::vector<std::size_t>{rows} : std::vector<std::size_t>{rows, width});
    parallel_for(rows, [&](std::size_t r) {
        double buf[DATA_CALC_MAX_COMPONENTS];
        const std::size_t w = calc_eval(rPool, root, rOperands, r, buf);
        for (std::size_t k = 0; k < width; ++k)
            detail::write_double(dst, r * width + k, buf[w == 1 ? 0 : k]);
    });
    return dst;
}

}  // namespace

Mesh data_calc(const Mesh& rMesh, const std::string& rExpression, const DataCalcOptions& rOpts) {
    if (rOpts.output.empty())
        calc_fail("an output array name is required");
    if (!rOpts.overwrite && data_has(rMesh, rOpts.location, rOpts.output))
        calc_fail("output name '" + rOpts.output + "' already exists in " +
                  data_location_name(rOpts.location) + " (pass overwrite=true to replace it)");

    const std::vector<CalcToken> tokens = calc_tokenize(rExpression);
    if (tokens.size() == 1)  // just the End token
        calc_fail("the expression is empty");

    std::vector<CalcNode> pool;
    CalcParser parser(tokens, pool);
    const std::size_t root = parser.ParseAll();

    std::vector<CalcOperand> operands;
    calc_infer(pool, root, rMesh, rOpts.location, operands);

    const std::size_t width = pool[root].mWidth;
    if (width > DATA_CALC_MAX_COMPONENTS)
        calc_fail("result has " + std::to_string(width) +
                  " components, more than the supported maximum of " +
                  std::to_string(DATA_CALC_MAX_COMPONENTS));

    Mesh out = detail::clone_mesh(rMesh);

    if (rOpts.location == DataLocation::Point) {
        const std::size_t rows = rMesh.NumPoints();
        calc_load_operands(
            operands, rows,
            [&](const std::string& n) -> const NDArray& { return rMesh.PointData(n); },
            "point_data");
        out.AddPointData(rOpts.output, calc_run(pool, root, operands, rows, width, rOpts.dtype));
        return out;
    }

    if (rOpts.location == DataLocation::Field) {
        // Field data has no entity correspondence, so the "row count" is simply
        // the common length of every referenced array.
        std::size_t rows = 0;
        for (const CalcOperand& op : operands) {
            const std::size_t r = detail::rows(rMesh.FieldData(op.mName));
            if (rows == 0)
                rows = r;
            else if (r != rows)
                calc_fail("field_data arrays disagree in length (" + std::to_string(rows) + " vs " +
                          std::to_string(r) + ")");
        }
        calc_load_operands(
            operands, rows,
            [&](const std::string& n) -> const NDArray& { return rMesh.FieldData(n); },
            "field_data");
        out.AddFieldData(rOpts.output, calc_run(pool, root, operands, rows, width, rOpts.dtype));
        return out;
    }

    // Cell location: evaluate once per block, producing exactly one array per
    // cell block as the uniform API requires.
    const std::size_t nblocks = rMesh.NumCellBlocks();
    for (const CalcOperand& op : operands)
        if (rMesh.CellDataNumBlocks(op.mName) != nblocks)
            calc_fail("cell_data '" + op.mName + "' has " +
                      std::to_string(rMesh.CellDataNumBlocks(op.mName)) +
                      " block(s) but the mesh has " + std::to_string(nblocks) + " cell block(s)");

    std::vector<NDArray> blocks;
    blocks.reserve(nblocks);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const std::size_t rows = rMesh.Cells(b).NumCells();
        calc_load_operands(
            operands, rows,
            [&](const std::string& n) -> const NDArray& { return rMesh.CellData(n, b); },
            "cell_data");
        blocks.push_back(calc_run(pool, root, operands, rows, width, rOpts.dtype));
    }
    out.AddCellData(rOpts.output, std::move(blocks));
    return out;
}

}  // namespace meshioplusplus
