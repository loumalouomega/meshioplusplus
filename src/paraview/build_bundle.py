"""Build the self-contained ParaView plugin bundle for one ParaView Python.

    python src/paraview/build_bundle.py --platform linux-x64 --python-tag 312 \
        --local-wheel 'dist/meshioplusplus-*cp312*.whl' --outdir dist

The plugin's hard step is installing meshio++ into ParaView's own Python, and
Kitware's binaries make it harder (no pip, a Python of their own). The bundle
skips it: a zip holding the plugin beside a ``site-packages/`` directory with
meshio++ and its runtime dependencies already unpacked, which the plugin puts
on ``sys.path`` when it finds it. Unzip it anywhere and load the plugin file
(Tools > Manage Plugins > Load New, or ``PV_PLUGIN_PATH``); nothing is
installed. See doc/paraview_plugin.md.

Why wheels unpacked, not wheels. ``pvpython`` has no installer to hand them
to, and a wheel is a zip whose ``.py`` and compiled extension can be imported
only once extracted.

Why numpy is **not** bundled. ParaView ships it, and VTK's numpy support is
built against that copy; a second numpy ahead of it on ``sys.path`` is the
quickest way to break ``vtkmodules.numpy_interface``. The compiled core uses
numpy's C API only through pybind11 at run time, so ParaView's numpy serves.

Why ``rich`` is. ``meshioplusplus/__init__.py`` imports ``_cli`` eagerly, and
``_cli`` imports rich (the Blender extension's reasoning, src/blender).

The add-ons (``meshioplusplus-occt``, ``-netgen``, ``-mmg``) go in when their
wheels are given (``--addon-wheel``) or can be downloaded for the platform;
an add-on with no wheel for the platform is left out with a message.
"""

from __future__ import annotations

import argparse
import glob
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import zipfile

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
PLUGIN = REPO / "tools" / "paraview-meshioplusplus-plugin.py"

#: Platform -> the PEP 425 platform tags to ask pip for, most specific first
#: (the Blender extension's table).
PLATFORMS = {
    "linux-x64": ["manylinux_2_34_x86_64", "manylinux2014_x86_64"],
    "windows-x64": ["win_amd64"],
    "macos-x64": ["macosx_13_0_x86_64", "macosx_11_0_x86_64", "macosx_10_9_x86_64"],
    "macos-arm64": ["macosx_14_0_arm64", "macosx_13_0_arm64", "macosx_11_0_arm64"],
}
#: Pure-Python runtime requirements; pip resolves rich's own.
PURE_REQUIREMENTS = ["rich"]
ADDONS = ("meshioplusplus-occt", "meshioplusplus-netgen", "meshioplusplus-mmg")
#: Distributions never bundled: ParaView ships numpy.
EXCLUDED = ("numpy",)


def project_version() -> str:
    with (REPO / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def _pip_download(dest, requirement, python_tag, platform_tags, deps=True):
    command = [
        sys.executable,
        "-m",
        "pip",
        "download",
        "--only-binary=:all:",
        "--dest",
        str(dest),
        "--python-version",
        python_tag,
        "--implementation",
        "cp",
    ]
    if not deps:
        command.append("--no-deps")
    for tag in platform_tags:
        command += ["--platform", tag]
    command.append(requirement)
    return subprocess.run(command, check=False).returncode == 0


def _copy_globs(patterns, dest, required=True):
    found = []
    for pattern in patterns or []:
        matches = sorted(glob.glob(pattern))
        if not matches and required:
            raise SystemExit(f"matched nothing: {pattern}")
        if not matches:
            print(f"build_bundle: {pattern} matched nothing; left out")
        for match in matches:
            shutil.copy2(match, dest / pathlib.Path(match).name)
            found.append(pathlib.Path(match).name)
    return found


def collect_wheels(dest, version, platform, python_tag, local_wheel, addon_wheels):
    """Every wheel the bundle unpacks, in ``dest``."""
    dest.mkdir(parents=True, exist_ok=True)
    tags = PLATFORMS[platform]
    if local_wheel:
        _copy_globs([local_wheel], dest)
    elif not _pip_download(dest, f"meshioplusplus=={version}", python_tag, tags, False):
        raise SystemExit(
            f"no meshioplusplus {version} wheel for cp{python_tag} {platform}"
        )
    given = _copy_globs(addon_wheels, dest, required=False)
    for addon in ADDONS:
        prefix = addon.replace("-", "_")
        if any(name.startswith(prefix) for name in given):
            continue
        if not _pip_download(dest, f"{addon}=={version}", python_tag, tags, False):
            print(f"build_bundle: no {addon} wheel for {platform}; left out")
    for requirement in PURE_REQUIREMENTS:
        if not _pip_download(dest, requirement, python_tag, ["any"]):
            raise SystemExit(f"could not download {requirement}")
    wheels = sorted(dest.glob("*.whl"))
    names = [w.name for w in wheels]
    if not any(n.startswith("meshioplusplus-") for n in names):
        raise SystemExit(f"no meshioplusplus wheel in {dest}")
    if not any(n.startswith("rich-") for n in names):
        raise SystemExit(
            "rich must be bundled: meshioplusplus/__init__.py imports _cli"
        )
    return [w for w in wheels if not w.name.startswith(EXCLUDED)]


def unpack_wheel(wheel, site_packages):
    """Unpack ``wheel`` as an installer would: its ``.data/purelib`` and
    ``.data/platlib`` into ``site_packages`` with the rest; its scripts,
    headers and data (the plugin file among them) left out."""
    with zipfile.ZipFile(wheel) as archive:
        for info in archive.infolist():
            name = info.filename
            parts = name.split("/")
            if parts[0].endswith(".data"):
                if len(parts) < 3 or parts[1] not in ("purelib", "platlib"):
                    continue
                name = "/".join(parts[2:])
            if not name or name.endswith("/"):
                continue
            target = site_packages / name
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)
            mode = (info.external_attr >> 16) & 0o777
            if mode:
                os.chmod(target, mode)


_README = """meshio++ {version} ParaView plugin, for ParaView's Python 3.{minor} ({platform})

1. Unzip this folder anywhere and keep it whole: the plugin finds meshio++ in
   the site-packages folder beside it.
2. In ParaView: Tools > Manage Plugins > Load New, choose
   paraview-meshioplusplus-plugin.py, and tick Auto Load. Or set
   PV_PLUGIN_PATH to this folder before starting ParaView.
3. From pvpython: paraview.simple.LoadPlugin(path, ns=globals()).

It needs a ParaView whose Python is 3.{minor} (pvpython -c "import sys;
print(sys.version)"). See https://loumalouomega.github.io/meshioplusplus/paraview_plugin
"""


def build(platform, outdir, python_tag, local_wheel=None, addon_wheels=None):
    version = project_version()
    root = f"meshioplusplus-paraview-{version}-{platform}-cp{python_tag}"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        wheels = collect_wheels(
            tmp / "wheels", version, platform, python_tag, local_wheel, addon_wheels
        )
        site = tmp / root / "site-packages"
        site.mkdir(parents=True)
        for wheel in wheels:
            unpack_wheel(wheel, site)
        shutil.copy2(PLUGIN, tmp / root / PLUGIN.name)
        (tmp / root / "README.txt").write_text(
            _README.format(version=version, minor=python_tag[1:], platform=platform)
        )
        for name in ("LICENSE", "LICENSE-BINARIES.md", "THIRD_PARTY_NOTICES.md"):
            if (REPO / name).exists():
                shutil.copy2(REPO / name, tmp / root / name)

        outdir.mkdir(parents=True, exist_ok=True)
        target = outdir / f"{root}.zip"
        # Deterministic: sorted names, a fixed timestamp, the files' modes.
        fixed = (1980, 1, 1, 0, 0, 0)
        with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted((tmp / root).rglob("*")):
                if path.is_dir():
                    continue
                info = zipfile.ZipInfo(
                    f"{root}/{path.relative_to(tmp / root).as_posix()}", date_time=fixed
                )
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = (path.stat().st_mode & 0o777 or 0o644) << 16
                archive.writestr(info, path.read_bytes())
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", choices=sorted(PLATFORMS), required=True)
    parser.add_argument(
        "--python-tag",
        required=True,
        help="the CPython of the target ParaView, e.g. 312 for Python 3.12",
    )
    parser.add_argument(
        "--local-wheel",
        help="glob for a locally built meshioplusplus wheel (CI, or a version "
        "not on PyPI yet); without it the release wheel is downloaded",
    )
    parser.add_argument(
        "--addon-wheel",
        action="append",
        help="glob for an add-on wheel to bundle (repeatable; one that "
        "matches nothing is left out); add-ons not given are downloaded when "
        "they exist for the platform",
    )
    parser.add_argument("--outdir", type=pathlib.Path, default=HERE / "dist")
    args = parser.parse_args(argv)
    target = build(
        args.platform, args.outdir, args.python_tag, args.local_wheel, args.addon_wheel
    )
    print(f"wrote {target} ({os.path.getsize(target) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
