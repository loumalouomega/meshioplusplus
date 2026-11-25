"""Verify CycloneDX SBOM generation and THIRD_PARTY_NOTICES.md pin parity (roadmap 3.14.6)."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[2]
_SPEC = importlib.util.spec_from_file_location(
    "gen_sbom", REPO / "tools" / "gen_sbom.py"
)
assert _SPEC and _SPEC.loader
sbom = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(sbom)


def _notices_table() -> dict[str, str]:
    """Parse ``| Component | Version | ... |`` rows from ``THIRD_PARTY_NOTICES.md``."""
    text = (REPO / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    versions: dict[str, str] = {}
    for line in text.splitlines():
        if (
            not line.startswith("|")
            or line.startswith("| Component")
            or line.startswith("| ---")
        ):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2:
            continue
        # Strip markdown link syntax `[Label](url)` to plain name.
        name = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", cells[0])
        versions[name] = cells[1]
    return versions


def test_third_party_notices_matches_pinned_versions():
    notices = _notices_table()
    pins = sbom.collect_pins(REPO)
    # Map pin key -> substring identifying the component row in THIRD_PARTY_NOTICES.md.
    expected_rows = {
        "opencascade": "Open CASCADE Technology",
        "zlib": "zlib",
        "hdf5": "HDF5",
        "netcdf-c": "netCDF-C",
        "cgns": "CGNS",
        "zstd": "zstd",
        "lz4": "LZ4",
        "bzip2": "bzip2",
        "kahip": "KaHIP",
        "cutcells": "CutCells",
        "netgen": "Netgen",
        "mmg": "MMG",
        "pugixml": "pugixml",
        "gidpost": "gidpost",
        "nlohmann-json": "nlohmann/json",
    }
    for key, label in expected_rows.items():
        matches = [ver for row_name, ver in notices.items() if label in row_name]
        assert matches, f"component {label!r} missing from THIRD_PARTY_NOTICES.md"
        pinned = pins[key]["version"]
        assert matches[0].startswith(pinned), (
            f"THIRD_PARTY_NOTICES.md has {label} = {matches[0]!r}, "
            f"but pin scripts specify {pinned!r}"
        )


def test_sbom_structure_and_artifact_filtering(tmp_path):
    wheel_bom = sbom.build_sbom(artifact="wheel", repo=REPO)
    assert wheel_bom["bomFormat"] == "CycloneDX"
    assert wheel_bom["specVersion"] == "1.5"
    wheel_names = {c["name"] for c in wheel_bom["components"]}
    assert {
        "HDF5",
        "netCDF-C",
        "CGNS",
        "zlib",
        "zstd",
        "LZ4",
        "KaHIP",
        "CutCells",
    } <= wheel_names
    assert "Open CASCADE Technology" not in wheel_names

    cli_bom = sbom.build_sbom(artifact="cli", repo=REPO)
    cli_names = {c["name"] for c in cli_bom["components"]}
    assert "Open CASCADE Technology" in cli_names
    assert "HDF5" not in cli_names

    # CMakeCache filtering drops disabled components.
    cache = tmp_path / "CMakeCache.txt"
    cache.write_text(
        "MESHIOPLUSPLUS_WITH_HDF5:BOOL=OFF\nMESHIOPLUSPLUS_WITH_ZLIB:BOOL=ON\n"
    )
    filtered = sbom.build_sbom(artifact="wheel", repo=REPO, cmake_cache=cache)
    filtered_names = {c["name"] for c in filtered["components"]}
    assert "HDF5" not in filtered_names
    assert "zlib" in filtered_names

    # CLI entrypoint writes valid JSON.
    out = tmp_path / "bom.cdx.json"
    assert sbom.main(["--artifact", "wheel", "--out", str(out)]) == 0
    parsed = json.loads(out.read_text(encoding="utf-8"))
    assert parsed["bomFormat"] == "CycloneDX"
