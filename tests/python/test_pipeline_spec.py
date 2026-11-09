"""The pipeline step spec holds both engines (roadmap §7).

``tools/pipeline_steps.json`` lists every step, its keys and, for the keys both
engines read with a literal default, that key's kind (``number``, ``flag``,
``text``) and default. ``tools/gen_pipeline_table.py`` writes the two step
tables from it. These tests check that the generated tables are current and
that every literal-default read of a step's option, in ``_pipeline.py``'s
``_apply_step`` and in ``pipeline.cpp``'s ``pipe_apply_step``, names a key of
that step with the spec's kind and default. A default changed in one engine
only, or a key read under the wrong kind, fails here.
"""

import ast
import json
import pathlib
import re
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SPEC = ROOT / "tools" / "pipeline_steps.json"
PY_SRC = ROOT / "src" / "python" / "meshioplusplus" / "_pipeline.py"
CPP_SRC = ROOT / "src" / "cpp" / "src" / "operations" / "pipeline.cpp"

pytestmark = pytest.mark.skipif(
    not (SPEC.exists() and PY_SRC.exists() and CPP_SRC.exists()),
    reason="needs a source checkout",
)

_PY_KINDS = {"_number": "number", "_flag": "flag", "_text": "text"}


def _spec_keys():
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    out = {}
    for st in spec["steps"] + spec["v2_steps"]:
        out[st["op"]] = {k["name"]: k for k in st["keys"]}
    return out


def _same(kind, default, value):
    if kind == "number":
        return isinstance(value, (int, float)) and float(value) == float(default)
    return type(value) is type(default) and value == default


def _python_reads():
    """(ops, key, kind, default, line) per literal-default read in _apply_step."""
    tree = ast.parse(PY_SRC.read_text(encoding="utf-8"))
    fn = next(
        f
        for f in ast.walk(tree)
        if isinstance(f, ast.FunctionDef) and f.name == "_apply_step"
    )

    def ops_of(test):
        if not (
            isinstance(test, ast.Compare)
            and isinstance(test.left, ast.Name)
            and test.left.id == "op"
        ):
            return []
        c = test.comparators[0]
        if isinstance(test.ops[0], ast.Eq) and isinstance(c, ast.Constant):
            return [c.value]
        if isinstance(test.ops[0], ast.In) and isinstance(c, (ast.Tuple, ast.List)):
            return [e.value for e in c.elts if isinstance(e, ast.Constant)]
        return []

    reads = []
    for node in ast.walk(fn):
        if not isinstance(node, ast.If):
            continue
        ops = tuple(ops_of(node.test))
        if ops:
            for st in node.body:
                for c in ast.walk(st):
                    if not (
                        isinstance(c, ast.Call)
                        and isinstance(c.func, ast.Name)
                        and c.func.id in _PY_KINDS
                        and len(c.args) >= 3
                        and isinstance(c.args[0], ast.Name)
                        and c.args[0].id == "step"
                        and isinstance(c.args[1], ast.Constant)
                    ):
                        continue
                    try:
                        default = ast.literal_eval(c.args[2])
                    except ValueError:
                        continue  # a computed default: nothing to pin
                    reads.append(
                        (ops, c.args[1].value, _PY_KINDS[c.func.id], default, c.lineno)
                    )
    return reads


def _cpp_reads():
    """(ops, key, kind, default, line) per literal-default read in pipe_apply_step."""
    text = CPP_SRC.read_text(encoding="utf-8")
    reads = []
    for m in re.finditer(r'if \(op == "\w+"(?: \|\| op == "\w+")*\) \{', text):
        ops = tuple(re.findall(r'op == "(\w+)"', m.group(0)))
        i, depth = m.end(), 1
        while depth:
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            i += 1
        body = text[m.end() : i]
        for r in re.finditer(
            r'pipe_(number|flag|text)\(rStep, "(\w+)", ([^()]*?)\)', body
        ):
            kind, key, raw = r.groups()
            raw = raw.strip()
            if raw in ("true", "false"):
                default = raw == "true"
            else:
                try:
                    default = ast.literal_eval(raw)
                except (ValueError, SyntaxError):
                    continue  # a struct member's default: nothing to pin
            line = text.count("\n", 0, m.end() + r.start()) + 1
            reads.append((ops, key, kind, default, line))
    return reads


def test_the_generated_tables_are_current():
    result = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "gen_pipeline_table.py"), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("engine", ["python", "cpp"])
def test_every_literal_default_read_matches_the_spec(engine):
    spec = _spec_keys()
    reads = _python_reads() if engine == "python" else _cpp_reads()
    assert len(reads) > 300, len(reads)  # the scan still sees the steps
    wrong = []
    for ops, key, kind, default, line in reads:
        # A block several steps share reads keys of one of them only.
        owners = [op for op in ops if key in spec.get(op, {})]
        if not owners:
            wrong.append(
                f"{'/'.join(ops)}.{key} (line {line}) is not a key of the step"
            )
            continue
        op = owners[0]
        entry = spec[op][key]
        if "kind" not in entry:
            continue  # read with several defaults by design (Generate's shapes)
        if engine == "cpp" and kind == "text" and default == entry.get("cpp_default"):
            continue  # C++ passes "" and the name lookup picks the same default
        if kind != entry["kind"] or not _same(kind, entry["default"], default):
            wrong.append(
                f"{op}.{key} (line {line}): {kind} {default!r}, "
                f"spec {entry['kind']} {entry['default']!r}"
            )
    assert not wrong, "\n".join(wrong)


def test_every_spec_default_is_read_by_an_engine():
    read = {(op, key) for ops, key, *_ in _python_reads() + _cpp_reads() for op in ops}
    stale = [
        f"{op}.{key}"
        for op, keys in _spec_keys().items()
        for key, entry in keys.items()
        if "kind" in entry and (op, key) not in read
    ]
    assert not stale, stale
