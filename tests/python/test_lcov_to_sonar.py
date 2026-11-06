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
    inside = root / "src" / "cpp" / "src" / "a.cpp"
    return f"""TN:
SF:{inside}
DA:1,3
DA:2,0
DA:4,1
BRDA:4,0,0,2
BRDA:4,0,1,-
end_of_record
SF:{root.parent / "elsewhere" / "gen.cpp"}
DA:1,1
end_of_record
SF:{inside}
DA:2,5
BRDA:4,0,1,1
end_of_record
"""


def test_converts_merges_and_drops_outside_files(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    trace = tmp_path / "cov.info"
    trace.write_text(_tracefile(root))
    out = tmp_path / "sonar.xml"
    assert lcov_to_sonar.main([str(trace), str(out), "--root", str(root)]) == 0

    doc = ET.parse(out).getroot()
    assert doc.tag == "coverage" and doc.get("version") == "1"
    files = doc.findall("file")
    assert [f.get("path") for f in files] == ["src/cpp/src/a.cpp"]
    lines = {
        int(e.get("lineNumber")): e.attrib for e in files[0].findall("lineToCover")
    }
    assert sorted(lines) == [1, 2, 4]
    assert lines[1]["covered"] == "true"
    # Line 2 was missed by the first record and hit by the second: merged hits.
    assert lines[2]["covered"] == "true"
    assert "branchesToCover" not in lines[1]
    # Branch (0,1) was "-" in one record and taken in the other.
    assert lines[4]["branchesToCover"] == "2"
    assert lines[4]["coveredBranches"] == "2"


def test_uncovered_lines_and_branches(tmp_path):
    files = lcov_to_sonar.parse_lcov(
        ["SF:src/x.cpp", "DA:7,0", "BRDA:7,0,0,-", "BRDA:7,0,1,0", "end_of_record"]
    )
    doc = ET.fromstring(lcov_to_sonar.to_sonar_xml(files, str(tmp_path)))
    (line,) = doc.find("file").findall("lineToCover")
    assert line.attrib == {
        "lineNumber": "7",
        "covered": "false",
        "branchesToCover": "2",
        "coveredBranches": "0",
    }
