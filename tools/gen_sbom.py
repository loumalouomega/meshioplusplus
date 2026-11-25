#!/usr/bin/env python3
"""Generate a CycloneDX 1.5 SBOM for meshio++ artifacts (roadmap 3.14.6).

Reads the pinned upstream library versions from ``build/deps-pins.sh``,
``build/occt-recipe.sh`` and ``CMakeLists.txt`` (zlib), plus the meshio++
release version from ``pyproject.toml``, and optionally filters/enriches
components from a build's ``CMakeCache.txt``.

Usage::

    python3 tools/gen_sbom.py --artifact wheel --out wheel.cdx.json
    python3 tools/gen_sbom.py --artifact cli --out cli.cdx.json
    python3 tools/gen_sbom.py --artifact all --cmake-cache build/CMakeCache.txt
"""

from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]

_ASSIGN_RE = re.compile(r'^\s*([A-Z0-9_]+)="([^"]+)"', re.M)
_ZLIB_URL_RE = re.compile(r"zlib-(\d+\.\d+\.\d+)\.tar\.gz")
_ZLIB_HASH_RE = re.compile(
    r"meshioplusplus_zlib[\s\S]*?URL_HASH\s+SHA256=([0-9a-f]{64})"
)
_PYPROJECT_VER_RE = re.compile(r'^version\s*=\s*"([^"]+)"', re.M)


def parse_shell_pins(path: pathlib.Path) -> dict[str, str]:
    """Extract ``KEY="value"`` assignments from a shell pin file."""
    return dict(_ASSIGN_RE.findall(path.read_text(encoding="utf-8")))


def parse_occt_pin(path: pathlib.Path, pin: str = "7.6") -> dict[str, str]:
    """Extract the OCCT version, tag, and commit for ``pin`` (7.6 or 7.9)."""
    text = path.read_text(encoding="utf-8")
    block = re.search(rf"^\s*{re.escape(pin)}\)\s*$([\s\S]*?)^\s*;;", text, re.M)
    if not block:
        raise ValueError(f"OCCT pin {pin!r} not found in {path}")
    return dict(_ASSIGN_RE.findall(block.group(1)))


def parse_zlib_pin(cmake_lists: pathlib.Path) -> tuple[str, str]:
    """Extract the static zlib version and SHA-256 from ``CMakeLists.txt``."""
    text = cmake_lists.read_text(encoding="utf-8")
    ver = _ZLIB_URL_RE.search(text)
    sha = _ZLIB_HASH_RE.search(text)
    if not ver or not sha:
        raise ValueError(f"zlib static pin not found in {cmake_lists}")
    return ver.group(1), sha.group(1)


def parse_project_version(pyproject: pathlib.Path) -> str:
    match = _PYPROJECT_VER_RE.search(pyproject.read_text(encoding="utf-8"))
    if not match:
        raise ValueError(f"project version not found in {pyproject}")
    return match.group(1)


def parse_cmake_cache(cache_path: pathlib.Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in cache_path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"([^:#/]+):[^=]+=(.*)", line)
        if match:
            values[match.group(1)] = match.group(2)
    return values


def collect_pins(
    repo: pathlib.Path = REPO, occt_pin: str = "7.6"
) -> dict[str, dict[str, str]]:
    """Return every pinned dependency with version, licence, sha256/commit, and scope."""
    deps = parse_shell_pins(repo / "build" / "deps-pins.sh")
    occt = parse_occt_pin(repo / "build" / "occt-recipe.sh", pin=occt_pin)
    zlib_ver, zlib_sha = parse_zlib_pin(repo / "CMakeLists.txt")

    return {
        "opencascade": {
            "name": "Open CASCADE Technology",
            "version": occt["OCCT_VERSION"],
            "license": "LGPL-2.1-only WITH OCCT-exception-1.0",
            "commit": occt["OCCT_COMMIT"],
            "url": "https://github.com/Open-Cascade-SAS/OCCT",
            "artifacts": "cli,occt-wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_OCCT",
        },
        "zlib": {
            "name": "zlib",
            "version": zlib_ver,
            "license": "Zlib",
            "sha256": zlib_sha,
            "url": "https://zlib.net/",
            "artifacts": "wheel,cli,wasm",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_ZLIB",
        },
        "hdf5": {
            "name": "HDF5",
            "version": deps["HDF5_VERSION"],
            "license": "BSD-3-Clause",
            "sha256": deps["HDF5_SHA256"],
            "commit": deps["HDF5_COMMIT"],
            "url": "https://www.hdfgroup.org/solutions/hdf5/",
            "artifacts": "wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_HDF5",
        },
        "netcdf-c": {
            "name": "netCDF-C",
            "version": deps["NETCDF_VERSION"],
            "license": "BSD-3-Clause",
            "sha256": deps["NETCDF_SHA256"],
            "commit": deps["NETCDF_COMMIT"],
            "url": "https://www.unidata.ucar.edu/software/netcdf/",
            "artifacts": "wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_NETCDF",
        },
        "cgns": {
            "name": "CGNS",
            "version": deps["CGNS_VERSION"],
            "license": "Zlib",
            "sha256": deps["CGNS_SHA256"],
            "commit": deps["CGNS_COMMIT"],
            "url": "https://cgns.github.io/",
            "artifacts": "wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_CGNS",
        },
        "zstd": {
            "name": "zstd",
            "version": deps["ZSTD_VERSION"],
            "license": "BSD-3-Clause",
            "sha256": deps["ZSTD_SHA256"],
            "commit": deps["ZSTD_COMMIT"],
            "url": "https://facebook.github.io/zstd/",
            "artifacts": "wheel",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_ZSTD",
        },
        "lz4": {
            "name": "LZ4",
            "version": deps["LZ4_VERSION"],
            "license": "BSD-2-Clause",
            "sha256": deps["LZ4_SHA256"],
            "commit": deps["LZ4_COMMIT"],
            "url": "https://lz4.org/",
            "artifacts": "wheel",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_LZ4",
        },
        "bzip2": {
            "name": "bzip2",
            "version": "1.0.8",
            "license": "bzip2-1.0.6",
            "url": "https://sourceware.org/bzip2/",
            "artifacts": "wheel",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_BZIP2",
        },
        "kahip": {
            "name": "KaHIP",
            "version": deps["KAHIP_VERSION"],
            "license": "MIT",
            "sha256": deps["KAHIP_SHA256"],
            "commit": deps["KAHIP_COMMIT"],
            "url": "https://kahip.github.io/",
            "artifacts": "wheel",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_KAHIP",
        },
        "cutcells": {
            "name": "CutCells",
            "version": deps["CUTCELLS_VERSION"],
            "license": "MIT",
            "commit": deps["CUTCELLS_COMMIT"],
            "url": "https://github.com/sclaus2/CutCells",
            "artifacts": "wheel",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_CUTCELLS",
        },
        "netgen": {
            "name": "Netgen",
            "version": deps["NETGEN_VERSION"],
            "license": "LGPL-2.1-only",
            "commit": deps["NETGEN_COMMIT"],
            "url": "https://ngsolve.org/",
            "artifacts": "netgen-wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_NETGEN",
        },
        "mmg": {
            "name": "MMG",
            "version": deps["MMG_VERSION"],
            "license": "LGPL-3.0-or-later",
            "commit": deps["MMG_COMMIT"],
            "url": "https://www.mmgtools.org/",
            "artifacts": "mmg-wheel,wasm-full",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_MMG",
        },
        "pugixml": {
            "name": "pugixml",
            "version": "1.14",
            "license": "MIT",
            "url": "https://pugixml.org/",
            "artifacts": "wheel,cli,wasm",
        },
        "gidpost": {
            "name": "gidpost",
            "version": "2.14",
            "license": "BSD-2-Clause-Views",
            "url": "https://www.gidsimulation.com/",
            "artifacts": "wheel,cli,wasm",
            "cmake_flag": "MESHIOPLUSPLUS_WITH_GIDPOST",
        },
        "nlohmann-json": {
            "name": "nlohmann/json",
            "version": "3.12.0",
            "license": "MIT",
            "url": "https://github.com/nlohmann/json",
            "artifacts": "wheel,cli,wasm",
        },
    }


def build_sbom(
    artifact: str = "all",
    repo: pathlib.Path = REPO,
    occt_pin: str = "7.6",
    cmake_cache: pathlib.Path | None = None,
) -> dict:
    """Build a CycloneDX 1.5 BOM dictionary."""
    project_version = parse_project_version(repo / "pyproject.toml")
    pins = collect_pins(repo, occt_pin=occt_pin)
    cache = parse_cmake_cache(cmake_cache) if cmake_cache else None

    components = []
    for key, info in sorted(pins.items()):
        targets = info["artifacts"].split(",")
        if artifact != "all" and artifact not in targets:
            continue
        if cache is not None and "cmake_flag" in info:
            flag = cache.get(info["cmake_flag"], "")
            if flag.upper() in ("OFF", "FALSE", "0"):
                continue
        comp: dict = {
            "type": "library",
            "bom-ref": f"pkg:generic/{key}@{info['version']}",
            "name": info["name"],
            "version": info["version"],
            "licenses": [{"expression": info["license"]}],
            "externalReferences": [{"type": "website", "url": info["url"]}],
        }
        if "sha256" in info:
            comp["hashes"] = [{"alg": "SHA-256", "content": info["sha256"]}]
        if "commit" in info:
            comp["properties"] = [{"name": "git:commit", "value": info["commit"]}]
        components.append(comp)

    return {
        "bomFormat": "CycloneDX",
        "specVersion": "1.5",
        "version": 1,
        "metadata": {
            "component": {
                "type": "library" if artifact != "cli" else "application",
                "name": "meshioplusplus",
                "version": project_version,
                "licenses": [{"expression": "MIT"}],
                "properties": [{"name": "meshioplusplus:artifact", "value": artifact}],
            }
        },
        "components": components,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--artifact",
        choices=(
            "wheel",
            "cli",
            "wasm",
            "wasm-full",
            "occt-wheel",
            "netgen-wheel",
            "mmg-wheel",
            "all",
        ),
        default="all",
    )
    parser.add_argument("--occt-pin", default="7.6", choices=("7.6", "7.9"))
    parser.add_argument("--cmake-cache", type=pathlib.Path)
    parser.add_argument("--out", type=pathlib.Path)
    args = parser.parse_args(argv)

    bom = build_sbom(
        artifact=args.artifact,
        repo=REPO,
        occt_pin=args.occt_pin,
        cmake_cache=args.cmake_cache,
    )
    text = json.dumps(bom, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
        print(f"wrote {args.out} ({len(bom['components'])} components)")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
