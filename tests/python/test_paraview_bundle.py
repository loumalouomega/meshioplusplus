"""The self-contained ParaView bundle (src/paraview/build_bundle.py).

The builder is tested with stand-in wheels and no network. The bundle itself
is tested by ``test_bundle_loads_in_a_kitware_pvpython`` when CI hands it a
built bundle and a ParaView binary from Kitware (``MESHIOPLUSPLUS_PARAVIEW_BUNDLE``
and ``MESHIOPLUSPLUS_BUNDLE_PVPYTHON``): ``pvpython`` loads the plugin with
no meshio++ installed and nothing on ``PYTHONPATH``.
"""

import importlib.util
import json
import os
import pathlib
import subprocess
import zipfile

import numpy as np
import pytest

import meshioplusplus

REPO = pathlib.Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def builder():
    spec = importlib.util.spec_from_file_location(
        "build_bundle", REPO / "src" / "paraview" / "build_bundle.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _wheel(path, files):
    with zipfile.ZipFile(path, "w") as archive:
        for name, text in files.items():
            archive.writestr(name, text)
    return path


def test_unpack_wheel_merges_purelib_and_leaves_data_out(builder, tmp_path):
    wheel = _wheel(
        tmp_path / "pkg-1.0-py3-none-any.whl",
        {
            "pkg/__init__.py": "x = 1\n",
            "pkg-1.0.dist-info/METADATA": "Name: pkg\n",
            "pkg-1.0.data/purelib/extra/mod.py": "y = 2\n",
            "pkg-1.0.data/platlib/native/core.py": "z = 3\n",
            "pkg-1.0.data/data/share/pkg/plugin.py": "# data\n",
            "pkg-1.0.data/scripts/tool": "#!/bin/sh\n",
        },
    )
    site = tmp_path / "site"
    builder.unpack_wheel(wheel, site)
    found = sorted(
        p.relative_to(site).as_posix() for p in site.rglob("*") if p.is_file()
    )
    assert found == [
        "extra/mod.py",
        "native/core.py",
        "pkg-1.0.dist-info/METADATA",
        "pkg/__init__.py",
    ]


def test_build_lays_out_the_bundle(builder, tmp_path, monkeypatch):
    version = builder.project_version()
    local = _wheel(
        tmp_path / f"meshioplusplus-{version}-cp312-cp312-linux_x86_64.whl",
        {
            "meshioplusplus/__init__.py": "",
            f"meshioplusplus-{version}.dist-info/RECORD": "",
        },
    )

    def fake_download(dest, requirement, python_tag, platform_tags, deps=True):
        if requirement == "rich":
            _wheel(dest / "rich-13.0.0-py3-none-any.whl", {"rich/__init__.py": ""})
            _wheel(dest / "numpy-2.0.0-cp312-cp312-any.whl", {"numpy/__init__.py": ""})
            return True
        return False  # no add-on wheels for this platform

    monkeypatch.setattr(builder, "_pip_download", fake_download)
    first = builder.build("linux-x64", tmp_path / "a", "312", str(local))
    second = builder.build("linux-x64", tmp_path / "b", "312", str(local))
    assert first.read_bytes() == second.read_bytes()  # deterministic
    root = f"meshioplusplus-paraview-{version}-linux-x64-cp312"
    assert first.name == root + ".zip"
    names = set(zipfile.ZipFile(first).namelist())
    assert f"{root}/paraview-meshioplusplus-plugin.py" in names
    assert f"{root}/README.txt" in names
    assert f"{root}/LICENSE" in names
    assert f"{root}/site-packages/meshioplusplus/__init__.py" in names
    assert f"{root}/site-packages/rich/__init__.py" in names
    # ParaView ships numpy; a second one would shadow it.
    assert not any("/site-packages/numpy/" in n for n in names)


def test_the_plugin_puts_a_sibling_site_packages_first():
    text = (REPO / "tools" / "paraview-meshioplusplus-plugin.py").read_text()
    bootstrap = text.index('"site-packages")')
    assert bootstrap < text.index("import meshioplusplus  #")


_BUNDLE = os.environ.get("MESHIOPLUSPLUS_PARAVIEW_BUNDLE")
_PVPYTHON = os.environ.get("MESHIOPLUSPLUS_BUNDLE_PVPYTHON")

_PROBE = """
import json, sys
import paraview.simple as pvs
from paraview import servermanager
pvs.LoadPlugin(sys.argv[1], ns=globals())
import meshioplusplus
reader = meshioreader(FileName=sys.argv[2])
reader.UpdatePipeline()
grid = servermanager.Fetch(reader)
print(json.dumps({"cells": grid.GetNumberOfCells(),
                  "from": meshioplusplus.__file__}))
"""


@pytest.mark.skipif(
    not (_BUNDLE and _PVPYTHON), reason="no bundle and Kitware pvpython given"
)
def test_bundle_loads_in_a_kitware_pvpython(tmp_path):
    """roadmap 6.11: pvpython from a Kitware tarball loads the bundle without
    pip, and reads a file through it."""
    with zipfile.ZipFile(_BUNDLE) as archive:
        archive.extractall(tmp_path / "bundle")
    (plugin,) = (tmp_path / "bundle").glob("*/paraview-meshioplusplus-plugin.py")
    src = tmp_path / "in.vtu"
    meshioplusplus.write(
        src,
        meshioplusplus.Mesh(
            np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            [("tetra", [[0, 1, 2, 3]])],
        ),
    )
    script = tmp_path / "probe.py"
    script.write_text(_PROBE)
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    run = subprocess.run(
        [_PVPYTHON, str(script), str(plugin), str(src)],
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
    )
    assert run.returncode == 0, run.stderr[-3000:]
    seen = json.loads(run.stdout.strip().splitlines()[-1])
    assert seen["cells"] == 1
    assert str(tmp_path / "bundle") in seen["from"]
