"""
A JSON document model for the JSON format readers and writers (Kratos CAD JSON
today): the Python twin of ``src/cpp/src/detail/json_value.{hpp,cpp}``.

Parsing is strict RFC 8259 plus line and block comments; objects come back as
:class:`JsonObject`, a list of ``(key, value)`` pairs in file order (a duplicate
key keeps both), arrays as lists, integer lexemes that fit an int64 as ``int``
and every other number as ``float``. Writing reproduces the C++ bytes: the
``dump(4)`` layout of nlohmann/json (and so of Kratos's
``PrettyPrintJsonString``) or compact, integers and integral floats within 2^53
as integers, any other float as the shortest of ``%.15g``, ``%.16g`` and
``%.17g`` that reads back exactly.
"""

import json
import math

from ._exceptions import ReadError

_MAX_DEPTH = 256
_INT64_MIN = -(2**63)
_INT64_MAX = 2**63 - 1


class JsonObject(list):
    """An object's members, a list of ``(key, value)`` pairs in file order."""

    def find(self, key, default=None):
        for k, v in self:
            if k == key:
                return v
        return default

    def __contains__(self, key):
        return any(k == key for k, _ in self)


def _strip_comments(text):
    if "/" not in text:
        return text
    out = []
    i, n = 0, len(text)
    start = 0
    while i < n:
        c = text[i]
        if c == '"':
            i += 1
            while i < n and text[i] != '"':
                i += 2 if text[i] == "\\" else 1
            i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            out.append(text[start:i])
            while i < n and text[i] != "\n":
                i += 1
            start = i
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            out.append(text[start:i])
            close = text.find("*/", i + 2)
            if close < 0:
                raise ReadError("meshio++: {fmt}: malformed JSON: unterminated comment")
            out.append(" ")
            i = close + 2
            start = i
        else:
            i += 1
    out.append(text[start:])
    return "".join(out)


def _parse_int(s):
    v = int(s)
    if s.startswith("-") and v == 0:
        return -0.0
    if v < _INT64_MIN or v > _INT64_MAX:
        return float(s)
    return v


def _reject_constant(name):
    raise ValueError(f"{name} is not JSON")


def _depth(value, level=0):
    if level > _MAX_DEPTH:
        return level
    deepest = level
    if isinstance(value, JsonObject):
        for _, v in value:
            deepest = max(deepest, _depth(v, level + 1))
    elif isinstance(value, list):
        for v in value:
            deepest = max(deepest, _depth(v, level + 1))
    return deepest


def parse(data, fmt):
    """Parse a whole JSON text (``bytes`` or ``str``); ``ReadError`` naming ``fmt``."""
    try:
        text = (
            data.decode("utf-8-sig") if isinstance(data, (bytes, bytearray)) else data
        )
        text = _strip_comments(text)
        value = json.loads(
            text,
            object_pairs_hook=JsonObject,
            parse_int=_parse_int,
            parse_constant=_reject_constant,
        )
    except ReadError as exc:
        raise ReadError(str(exc).replace("{fmt}", fmt)) from None
    except (ValueError, UnicodeDecodeError, RecursionError) as exc:
        raise ReadError(f"meshio++: {fmt}: malformed JSON: {exc}") from None
    if _depth(value) > _MAX_DEPTH:
        raise ReadError(
            f"meshio++: {fmt}: malformed JSON: nesting deeper than {_MAX_DEPTH}"
        )
    return value


def format_number(value):
    """A float's JSON spelling, as ``detail::append_json_number`` writes it."""
    value = float(value)
    if not math.isfinite(value):
        return "null"
    if (
        value == math.floor(value)
        and abs(value) < 9007199254740992.0
        and not (value == 0.0 and math.copysign(1.0, value) < 0)
    ):
        return str(int(value))
    for fmt in ("%.15g", "%.16g", "%.17g"):
        s = fmt % value
        if float(s) == value:
            break
    return s


_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def _string(s):
    out = ['"']
    for ch in s:
        e = _ESCAPES.get(ch)
        if e is not None:
            out.append(e)
        elif ord(ch) < 0x20:
            out.append("\\u%04x" % ord(ch))
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def _write(out, value, indent, level):
    def newline(lv):
        if indent > 0:
            out.append("\n" + " " * (indent * lv))

    if value is None:
        out.append("null")
    elif value is True:
        out.append("true")
    elif value is False:
        out.append("false")
    elif isinstance(value, int):
        out.append(str(value))
    elif isinstance(value, float):
        out.append(format_number(value))
    elif isinstance(value, str):
        out.append(_string(value))
    elif isinstance(value, JsonObject):
        if not value:
            out.append("{}")
            return
        out.append("{")
        for i, (k, v) in enumerate(value):
            if i:
                out.append(",")
            newline(level + 1)
            out.append(_string(k))
            out.append(": " if indent > 0 else ":")
            _write(out, v, indent, level + 1)
        newline(level)
        out.append("}")
    elif isinstance(value, (list, tuple)):
        if not value:
            out.append("[]")
            return
        out.append("[")
        for i, v in enumerate(value):
            if i:
                out.append(",")
            newline(level + 1)
            _write(out, v, indent, level + 1)
        newline(level)
        out.append("]")
    else:
        raise TypeError(f"cannot write {type(value).__name__} as JSON")


def dumps(value, indent=4):
    """The text of ``value``, ``dump(4)``-style when ``indent > 0``, else compact."""
    out = []
    _write(out, value, indent, 0)
    return "".join(out)
