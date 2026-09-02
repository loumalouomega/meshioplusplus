"""Build the Blender extension zip: the add-on plus the wheels it needs.

    python src/blender/build_extension.py --platform linux-x64 \
        --local-wheel dist/meshioplusplus-*.whl --outdir dist

Why wheels at all. ``import meshioplusplus`` hard-requires the compiled
``_core``: 45 of its format subpackages import it unconditionally at module
scope, so there is no pure-Python path to fall back to and no amount of
vendoring source would work. Blender 4.2's extension system installs bundled
wheels itself, offline, which is exactly the mechanism this needs.

Why ``rich`` is in the bundle. ``meshioplusplus/__init__.py`` imports ``_cli``
eagerly, and ``_cli`` imports rich -- so rich (and its own dependencies) is a
hard runtime requirement, not an extra.

Why numpy is **not**. Blender ships it, and a second numpy ahead of Blender's
own on ``sys.path`` is a well-known way to break ``bpy``'s numpy interop. If a
future Blender stops shipping it, adding it to ``_PURE_REQUIREMENTS`` is the
escape hatch.
"""

from __future__ import annotations

import argparse
import glob
import os
import pathlib
import shutil
import subprocess
import sys
import zipfile

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent

#: Blender platform -> the PEP 425 platform tags to ask pip for, most specific
#: first. macOS wheels are built against a deployment target that moves between
#: releases, so several candidates are offered rather than one guess.
PLATFORMS = {
    "linux-x64": ["manylinux_2_34_x86_64", "manylinux2014_x86_64"],
    "windows-x64": ["win_amd64"],
    "macos-x64": ["macosx_13_0_x86_64", "macosx_11_0_x86_64", "macosx_10_9_x86_64"],
    "macos-arm64": ["macosx_14_0_arm64", "macosx_13_0_arm64", "macosx_11_0_arm64"],
}

#: Pure-Python runtime requirements. ``rich`` pulls markdown-it-py, mdurl and
#: pygments; pip resolves them, so they are not spelled out here and cannot
#: drift when rich changes its own dependencies.
_PURE_REQUIREMENTS = ["rich"]

#: The CPython the target Blender embeds. 4.2 through 4.5 are all 3.11.
DEFAULT_PYTHON_TAG = "311"


def project_version() -> str:
    with (REPO / "pyproject.toml").open("rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def _pip_download(dest, requirement, python_tag, platform_tags):
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
    for tag in platform_tags:
        command += ["--platform", tag]
    command.append(requirement)
    subprocess.run(command, check=True)


def collect_wheels(dest, version, platform, python_tag, local_wheel=None):
    """Put every wheel the extension bundles into ``dest``."""
    dest.mkdir(parents=True, exist_ok=True)
    platform_tags = PLATFORMS[platform]

    if local_wheel:
        matches = sorted(glob.glob(local_wheel))
        if not matches:
            raise SystemExit(f"--local-wheel matched nothing: {local_wheel}")
        for match in matches:
            shutil.copy2(match, dest / pathlib.Path(match).name)
    else:
        _pip_download(dest, f"meshioplusplus=={version}", python_tag, platform_tags)

    # The pure-Python dependencies are the same file on every platform.
    for requirement in _PURE_REQUIREMENTS:
        _pip_download(dest, requirement, python_tag, ["any"])

    wheels = sorted(p.name for p in dest.glob("*.whl"))
    if not any(w.startswith("meshioplusplus-") for w in wheels):
        raise SystemExit(f"no meshioplusplus wheel landed in {dest}")
    if not any(w.startswith("rich-") for w in wheels):
        raise SystemExit(
            "rich is missing: meshioplusplus/__init__.py imports _cli eagerly, "
            "so it is a hard runtime requirement"
        )
    return wheels


def render_manifest(version, platform, wheels):
    template = (HERE / "blender_manifest.toml.in").read_text()
    listed = "[\n" + "".join(f'  "./wheels/{w}",\n' for w in wheels) + "]"
    return (
        template.replace("@VERSION@", version)
        .replace("@PLATFORM@", platform)
        .replace("@WHEELS@", listed)
    )


def build(platform, outdir, python_tag, local_wheel=None):
    version = project_version()
    wheel_dir = HERE / "wheels"
    if wheel_dir.exists():
        shutil.rmtree(wheel_dir)
    wheels = collect_wheels(wheel_dir, version, platform, python_tag, local_wheel)
    manifest = render_manifest(version, platform, wheels)
    (HERE / "blender_manifest.toml").write_text(manifest)

    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / f"meshioplusplus-{version}-{platform}.zip"
    # Deterministic: sorted names and a fixed timestamp, so two builds of the
    # same inputs produce the same bytes.
    fixed = (1980, 1, 1, 0, 0, 0)
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:

        def add(path, arcname):
            info = zipfile.ZipInfo(arcname, date_time=fixed)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, path.read_bytes())

        add(HERE / "blender_manifest.toml", "blender_manifest.toml")
        for source in sorted((HERE / "addon").glob("*.py")):
            add(source, source.name)
        for name in wheels:
            add(wheel_dir / name, f"wheels/{name}")
        add(REPO / "LICENSE", "LICENSE")
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", choices=sorted(PLATFORMS), required=True)
    parser.add_argument(
        "--local-wheel",
        help=(
            "glob for a locally built meshioplusplus wheel. Use this in CI and "
            "on any branch that bumps the version: the wheel the manifest must "
            "name does not exist on PyPI yet."
        ),
    )
    parser.add_argument("--outdir", type=pathlib.Path, default=HERE / "dist")
    parser.add_argument(
        "--python-tag",
        default=DEFAULT_PYTHON_TAG,
        help="CPython the target Blender embeds (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    target = build(args.platform, args.outdir, args.python_tag, args.local_wheel)
    print(f"wrote {target} ({os.path.getsize(target) / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
