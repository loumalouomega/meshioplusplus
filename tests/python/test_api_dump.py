"""The public-API dump (``tools/api/python_api.json``) and its comparison rule."""

import importlib.util
import json
import pathlib

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, REPO / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


check = _load("check_api_changes")


def _fn(*params):
    return {"kind": "function", "signature": {"parameters": list(params)}}


def _p(name, default=None, kind="POSITIONAL_OR_KEYWORD"):
    row = {"name": name, "kind": kind}
    if default is not None:
        row["default"] = default
    return row


def _dump(**names):
    return {"schema_version": 1, "names": names}


def test_identical_dumps_have_no_changes():
    d = _dump(f=_fn(_p("a"), _p("b", "1")))
    assert check.compare(d, d) == []


def test_additions_are_not_breaking():
    old = _dump(f=_fn(_p("a")))
    new = _dump(f=_fn(_p("a"), _p("b", "None")), g=_fn())
    assert check.compare(old, new) == []


@pytest.mark.parametrize(
    "new, fragment",
    [
        (_dump(), "removed from __all__"),
        (_dump(f=_fn(_p("b", "1"))), "'a' removed or renamed"),
        (_dump(f=_fn(_p("a"), _p("b", "2"))), "default of 'b' changed"),
        (_dump(f=_fn(_p("a"), _p("b"))), "default of 'b' changed"),
        (_dump(f=_fn(_p("a"), _p("b", "1"), _p("c"))), "new required parameter 'c'"),
        (_dump(f={"kind": "constant"}), "changed from function to constant"),
    ],
)
def test_incompatible_changes_are_found(new, fragment):
    old = _dump(f=_fn(_p("a"), _p("b", "1")))
    found = check.compare(old, new)
    assert found and any(fragment in line for line in found), found


def test_unreadable_signatures_are_skipped():
    old = _dump(f={"kind": "function", "signature": None})
    new = _dump(f=_fn(_p("x")))
    assert check.compare(old, new) == []


def test_class_methods_are_compared():
    def cls(method_params):
        return {
            "kind": "class",
            "signature": None,
            "members": {
                "m": {"kind": "method", "signature": {"parameters": method_params}}
            },
        }

    old = _dump(C=cls([_p("self"), _p("x", "1")]))
    assert (
        check.compare(old, _dump(C=cls([_p("self"), _p("x", "1"), _p("y", "2")]))) == []
    )
    assert check.compare(old, _dump(C=cls([_p("self")])))
    gone = {"kind": "class", "signature": None, "members": {}}
    assert any("C.m: removed" in c for c in check.compare(old, _dump(C=gone)))


def test_the_committed_dump_is_fresh():
    try:
        import meshioplusplus  # noqa: F401
        from meshioplusplus import _core  # noqa: F401
    except ImportError:
        pytest.skip("the extension is not built")
    gen = _load("gen_api_dump")
    committed = json.loads(gen.OUTPUT.read_text(encoding="utf-8"))
    current = gen.build(allow_unavailable=True)
    unavailable = {n for n, e in current["names"].items() if e["kind"] == "unavailable"}
    # Optional extras may be missing here; compare everything that did resolve.
    for name in unavailable:
        current["names"].pop(name, None)
        committed["names"].pop(name, None)
    assert current == committed, "API changed: run python tools/gen_api_dump.py"
