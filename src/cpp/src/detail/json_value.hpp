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
 * @file detail/json_value.hpp
 * @brief A small JSON document model, parser and writer for the readers and
 * writers of JSON formats (Kratos CAD JSON today).
 *
 * The core's other JSON user, the pipeline parser, links the vendored
 * nlohmann/json behind `MESHIOPLUSPLUS_WITH_JSON`; a format reader cannot
 * depend on an optional submodule (the single header, WASM and every build
 * must read the file), so this one is dependency-free:
 *
 *  - **Parsing** is strict RFC 8259 plus line and block comments, which
 *    Kratos's own `Parameters` accept. Nesting is bounded, numbers go through
 *    `detail::parse_double` (locale-independent), `\u` escapes (surrogate
 *    pairs included) become UTF-8, and every malformed input throws
 *    `ReadError` naming the format and the byte offset.
 *  - **Objects keep their member order** (a vector of pairs, not a map), so a
 *    rewrite keeps the file's layout; a duplicate key keeps both, and `Find`
 *    returns the first.
 *  - **Writing** is the layout nlohmann's `dump(4)` (and so Kratos's
 *    `PrettyPrintJsonString`) produces, or compact. A number that was an
 *    integer lexeme, or is integral and within 2^53, is written as an integer;
 *    any other as the shortest of `%.15g`, `%.16g` and `%.17g` that reads back
 *    exactly. The Python twin (`meshioplusplus/_json_text.py`) writes the
 *    same bytes.
 *
 * Core-private (not installed).
 */

// System includes
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace meshioplusplus::detail {

/// An object member, defined below JsonValue (which it holds by value).
struct JsonMember;

/// One JSON value; a plain aggregate, its members public.
struct JsonValue {
    enum class Kind : std::uint8_t { Null, Bool, Number, String, Array, Object };

    Kind mKind = Kind::Null;
    bool mBool = false;
    /// The value of a number.
    double mNumber = 0.0;
    /// Whether the number was written as an integer that fits an int64.
    bool mIsInteger = false;
    /// Its value then.
    std::int64_t mInteger = 0;
    /// A string's UTF-8 bytes.
    std::string mString;
    /// An array's items.
    std::vector<JsonValue> mItems;
    /// An object's members, in file order.
    std::vector<JsonMember> mMembers;

    bool IsNull() const noexcept { return mKind == Kind::Null; }
    bool IsBool() const noexcept { return mKind == Kind::Bool; }
    bool IsNumber() const noexcept { return mKind == Kind::Number; }
    bool IsString() const noexcept { return mKind == Kind::String; }
    bool IsArray() const noexcept { return mKind == Kind::Array; }
    bool IsObject() const noexcept { return mKind == Kind::Object; }

    /// The first member named @p Key of an object, or null.
    const JsonValue* Find(std::string_view Key) const noexcept;

    static JsonValue Null();
    static JsonValue Boolean(bool Value);
    static JsonValue Number(double Value);
    static JsonValue Integer(std::int64_t Value);
    static JsonValue String(std::string Value);
    static JsonValue Array();
    static JsonValue Object();

    /// Append to an array.
    JsonValue& Push(JsonValue Value);
    /// Append a member to an object (no check for an existing key).
    JsonValue& Add(std::string Key, JsonValue Value);
};

/// One member of an object: its key and value. A named struct, not a
/// `std::pair`: clang with libstdc++ cannot complete a vector of pairs holding
/// the type being defined.
struct JsonMember {
    std::string mKey;
    JsonValue mValue;

    JsonMember(std::string Key, JsonValue Value)
        : mKey(std::move(Key)), mValue(std::move(Value)) {}
};

/// Whether two values are equal: same kind, numbers by value, objects member
/// by member in order.
bool json_equal(const JsonValue& rA, const JsonValue& rB) noexcept;

/**
 * @brief Parse a whole JSON text.
 * @param Text the bytes (a UTF-8 BOM is skipped).
 * @param pFormat the format name for error messages.
 * @throws ReadError on malformed input, nesting deeper than 256, or trailing
 * content.
 */
JsonValue parse_json(std::string_view Text, const char* pFormat);

/// The text of @p rValue, `dump(4)`-style when @p Indent > 0, else compact.
std::string write_json(const JsonValue& rValue, int Indent = 4);

/// Append @p Value's JSON spelling (see the file comment) to @p rOut.
void append_json_number(std::string& rOut, double Value);

}  // namespace meshioplusplus::detail
