"""Tests for tools/lcov_to_sonar.py (doc/static_analysis.md)."""

import importlib.util
import pathlib
import xml.etree.ElementTree as ET

REPO = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "lcov_to_sonar", REPO / "tools/lcov_to_sonar.py"
)
lcov_to_sonar = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lcov_to_sonar)


def _tracefile(root):
    return f"""TN:
SF:{root}/src/cpp/src/a.cpp
FN:3,f
DA:3,5
DA:4,0
DA:7,2
BRDA:7,0,0,1
BRDA:7,0,1,-
BRDA:7,0,2,0
end_of_record
SF:/usr/include/c++/13/vector
DA:10,4
end_of_record
SF:{root}/src/cpp/src/a.cpp
DA:4,3
end_of_record
SF:src/cpp/src/b.hpp
DA:1,0
end_of_record
"""


def test_lines_and_branches(tmp_path, capsys):
    root = tmp_path / "repo"
    info = tmp_path / "cov.info"
    info.write_text(_tracefile(root), encoding="utf-8")
    out = tmp_path / "cov.xml"
    assert lcov_to_sonar.main([str(info), str(out), "--root", str(root)]) == 0
    assert "2 of 3 source files" in capsys.readouterr().out
    doc = ET.parse(out).getroot()
    assert doc.tag == "coverage" and doc.get("version") == "1"
    files = {f.get("path"): f for f in doc.findall("file")}
    # A system header outside the root is left out; paths are root-relative.
    assert sorted(files) == ["src/cpp/src/a.cpp", "src/cpp/src/b.hpp"]
    lines = {int(e.get("lineNumber")): e.attrib for e in files["src/cpp/src/a.cpp"]}
    assert lines[3]["covered"] == "true"
    # Two records of one file add their hits (a TU linked into two binaries).
    assert lines[4]["covered"] == "true"
    assert lines[7]["branchesToCover"] == "3" and lines[7]["coveredBranches"] == "1"
    assert "branchesToCover" not in lines[3]
    assert files["src/cpp/src/b.hpp"][0].get("covered") == "false"


def test_empty_tracefile(tmp_path):
    info = tmp_path / "empty.info"
    info.write_text("", encoding="utf-8")
    out = tmp_path / "out.xml"
    lcov_to_sonar.main([str(info), str(out), "--root", str(tmp_path)])
    assert ET.parse(out).getroot().findall("file") == []
