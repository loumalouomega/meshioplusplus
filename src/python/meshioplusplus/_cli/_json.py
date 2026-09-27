"""The CLI's one JSON emitter.

Every verb's ``--json`` goes through :func:`emit_json`, so the output is always
valid JSON: NaN and infinities become ``null`` (``json.dumps`` would print a
bare ``NaN``, which no strict parser accepts), numpy scalars and arrays become
plain numbers and lists, and tuples become lists. The native CLI's
``json_out.hpp`` follows the same rules, and ``tests/python/test_cli_json.py``
holds the two to one shape.
"""

from __future__ import annotations

import json
import math

import numpy as np


def to_jsonable(obj):
    """``obj`` with every value JSON can carry, recursively; non-finite -> ``None``."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    return obj


def emit_json(obj) -> None:
    """Print ``obj`` as indented, strictly valid JSON."""
    print(json.dumps(to_jsonable(obj), indent=2, allow_nan=False, ensure_ascii=False))
