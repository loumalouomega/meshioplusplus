"""Tests for tools/dup_report.py (doc/static_analysis.md)."""

import importlib.util
import json
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "dup_report", REPO / "tools/dup_report.py"
)
dup_report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dup_report)

CPP_BODY = """
double {name}(const double* a, const double* b, int n) {{
    double s = 0.0;
    for (int i = 0; i < n; ++i) {{
        if (a[i] > b[i])
            s += a[i] * {k} - b[i];
        else
            s -= b[i] * a[i] + {k};
    }}
    return s / static_cast<double>(n > 0 ? n : 1);
}}
"""


def test_finds_a_renamed_cpp_clone(tmp_path, capsys):
    (tmp_path / "a.cpp").write_text("// a\n" + CPP_BODY.format(name="alpha", k=2))
    (tmp_path / "b.cpp").write_text(
        "#include <x>\nint z;\n" + CPP_BODY.format(name="beta", k=7)
    )
    (tmp_path / "c.cpp").write_text("int lonely() { return 1; }\n")
    out = tmp_path / "r.json"
    rc = dup_report.main(
        [
            "--lang",
            "cpp",
            "--root",
            str(tmp_path),
            "--paths",
            "a.cpp",
            "b.cpp",
            "c.cpp",
            "--tokens",
            "20",
            "--json",
            str(out),
        ]
    )
    assert rc == 0
    result = json.loads(out.read_text())["cpp"]
    assert result["files_scanned"] == 3
    names = {f["file"] for f in result["files"]}
    assert names == {"a.cpp", "b.cpp"}
    region = result["regions"][0]
    assert (region["a"], region["b"]) == ("a.cpp", "b.cpp")
    assert region["a_lines"] == [3, 12] and region["b_lines"] == [4, 13]
    assert "duplicated" in capsys.readouterr().out


def test_python_clone_and_unique_code(tmp_path):
    body = (
        "def {n}(xs, k):\n"
        "    total = 0\n"
        "    for x in xs:\n"
        "        if x > k:\n"
        "            total += x * 2 - k\n"
        "        else:\n"
        "            total -= k\n"
        "    return total / max(len(xs), 1)\n"
    )
    (tmp_path / "p.py").write_text(body.format(n="f"))
    (tmp_path / "q.py").write_text("# comment\n" + body.format(n="g"))
    (tmp_path / "r.py").write_text("VALUE = {'a': 1}\n")
    paths = [str(tmp_path / n) for n in ("p.py", "q.py", "r.py")]
    result = dup_report.analyse(paths, dup_report.py_tokens, 15, 40, str(tmp_path))
    assert {f["file"] for f in result["files"]} == {"p.py", "q.py"}
    assert result["regions"][0]["b_lines"][0] == 2
    assert 0.0 < result["density"] < 100.0


def _table(rows):
    return "".join(f'    {{"{k}", {v}}},\n' for k, v in rows)


def test_literal_tables_of_one_shape_do_not_match(tmp_path):
    # Roadmap 8.1: once strings and numbers are normalised, two literal tables
    # of the same shape look identical although no entry agrees.
    a = [(f"a{i}", i) for i in range(40)]
    b = [(f"b{i}", 100 + i) for i in range(40)]
    (tmp_path / "a.cpp").write_text("const T kA[] = {\n" + _table(a) + "};\n")
    (tmp_path / "b.cpp").write_text("const T kB[] = {\n" + _table(b) + "};\n")
    paths = [str(tmp_path / n) for n in ("a.cpp", "b.cpp")]
    result = dup_report.analyse(paths, dup_report.cpp_tokens, 40, 40, str(tmp_path))
    assert result["files"] == []
    assert result["regions"] == []


def test_python_literal_dicts_do_not_match(tmp_path):
    body = "".join(f'    "{{p}}{i}": ("{{p}}x{i}", {i}),\n' for i in range(30))
    (tmp_path / "p.py").write_text("A = {\n" + body.format(p="p") + "}\n")
    (tmp_path / "q.py").write_text("B = {\n" + body.format(p="q") + "}\n")
    paths = [str(tmp_path / n) for n in ("p.py", "q.py")]
    result = dup_report.analyse(paths, dup_report.py_tokens, 40, 40, str(tmp_path))
    assert result["files"] == []


def test_table_window_rule():
    assert dup_report.is_table_window(("{", "S", ",", "N", "}", ","))
    assert dup_report.is_table_window(("-", "N", ".", "N"))
    # One token in ten may be a name: a table's declaration beside its rows.
    assert dup_report.is_table_window(("I",) + ("S", ",") * 5)
    assert not dup_report.is_table_window(("{", "I", ",", "N", "}"))
    assert not dup_report.is_table_window(("return", "N", ";"))


def test_code_beside_a_table_still_matches(tmp_path):
    table = _table([(f"k{i}", i) for i in range(30)])
    (tmp_path / "a.cpp").write_text(
        "const T kA[] = {\n" + table + "};\n" + CPP_BODY.format(name="f", k=1)
    )
    (tmp_path / "b.cpp").write_text(CPP_BODY.format(name="g", k=3))
    paths = [str(tmp_path / n) for n in ("a.cpp", "b.cpp")]
    result = dup_report.analyse(paths, dup_report.cpp_tokens, 20, 40, str(tmp_path))
    assert {f["file"] for f in result["files"]} == {"a.cpp", "b.cpp"}


def test_r_binding_c_sources_are_scanned():
    assert "bindings/r/meshioplusplus/src" in dup_report.CPP_ROOTS
    assert ".c" in dup_report.LANGS["cpp"][1]
    assert "node_order.cpp" in dup_report.EXCLUDED_FILES
