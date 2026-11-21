"""Dump the public Python API of ``meshioplusplus`` to ``tools/api/python_api.json``.

Every name of ``meshioplusplus.__all__`` gets its kind and, where Python can read
one, its signature (parameters with their kinds and defaults); a class also lists
its public methods the same way. ``tools/check_api_changes.py`` compares two such
dumps, so a removed name, a removed parameter or a changed default is noticed and
has to be declared with a ``**Breaking:**`` line in ``CHANGELOG.md``::

    python tools/gen_api_dump.py            # rewrite the dump
    python tools/gen_api_dump.py --check    # exit 1 if it is stale

Generate it from an editable install with the extensions the CI ``test`` job
installs (``pip install -e ".[all]"``), so every lazily loaded name resolves; a
name that cannot be imported is written as ``"kind": "unavailable"`` and the run
fails unless ``--allow-unavailable`` is given.
"""

from __future__ import annotations

import argparse
import inspect
import json
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
OUTPUT = REPO / "tools" / "api" / "python_api.json"

_ADDRESS = re.compile(r" at 0x[0-9a-fA-F]+")


def _default(value):
    """A stable text for a default: primitives by repr, anything else by type."""
    if value is inspect.Parameter.empty:
        return None
    if value is None or isinstance(value, (bool, int, float, str, bytes)):
        return repr(value)
    if isinstance(value, (tuple, list, set, frozenset, dict)) and not value:
        return repr(value)
    if isinstance(value, (tuple, list)) and all(
        v is None or isinstance(v, (bool, int, float, str)) for v in value
    ):
        return repr(value)
    return "<" + type(value).__name__ + ">"


def _signature(obj):
    """``{"parameters": [...], "return": ...}`` or ``None`` when unreadable."""
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return None
    params = []
    for p in sig.parameters.values():
        row = {"name": p.name, "kind": p.kind.name}
        d = _default(p.default)
        if d is not None:
            row["default"] = _ADDRESS.sub("", d)
        params.append(row)
    return {"parameters": params}


def _kind(obj):
    if inspect.ismodule(obj):
        return "module"
    if inspect.isclass(obj):
        return "class"
    if callable(obj):
        return "function"
    return "constant"


def _class_members(cls):
    members = {}
    for name, value in sorted(vars(cls).items()):
        if name.startswith("_") and name not in ("__init__", "__call__"):
            continue
        func = value
        if isinstance(value, (staticmethod, classmethod)):
            func = value.__func__
        elif isinstance(value, property):
            members[name] = {"kind": "property"}
            continue
        if callable(func):
            members[name] = {"kind": "method", "signature": _signature(func)}
    return members


def build(allow_unavailable=False):
    import meshioplusplus

    names = sorted(set(meshioplusplus.__all__))
    api = {}
    missing = []
    for name in names:
        try:
            obj = getattr(meshioplusplus, name)
        except Exception as exc:  # noqa: BLE001 - any import failure is reported
            api[name] = {"kind": "unavailable", "error": type(exc).__name__}
            missing.append(name)
            continue
        entry = {"kind": _kind(obj)}
        if entry["kind"] in ("function", "class"):
            entry["signature"] = _signature(obj)
        if entry["kind"] == "class":
            entry["members"] = _class_members(obj)
        api[name] = entry
    if missing and not allow_unavailable:
        raise SystemExit(
            "cannot import: " + ", ".join(missing) + " (install the extras, or use "
            "--allow-unavailable)"
        )
    return {"schema_version": 1, "names": api}


def render(doc):
    return json.dumps(doc, indent=1, sort_keys=True) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check", action="store_true", help="exit 1 if the dump is stale"
    )
    parser.add_argument("--allow-unavailable", action="store_true")
    args = parser.parse_args(argv)
    text = render(build(args.allow_unavailable))
    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.is_file() else ""
        if current != text:
            print(f"{OUTPUT.relative_to(REPO)} is stale: run tools/gen_api_dump.py")
            return 1
        return 0
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(REPO)} ({len(json.loads(text)['names'])} names)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
