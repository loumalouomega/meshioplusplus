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
