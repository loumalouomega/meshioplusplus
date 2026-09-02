# meshio++ for Blender

A Blender 4.2+ **extension** that puts all 43 of meshio++'s formats behind
`File > Import` and `File > Export`. Blender reads STL, OBJ and PLY and
essentially nothing else from the FEA world; this closes that.

User documentation lives at [`doc/blender.md`](../../doc/blender.md). This file
is about building and hacking on it.

## Layout

| Path | What it is |
| --- | --- |
| `addon/` | exactly what ends up at the zip root — the Blender-facing shell |
| `blender_manifest.toml.in` | the manifest template; the real one is generated |
| `build_extension.py` | collects the wheels, renders the manifest, writes the zip |

**The conversion itself is not here.** It lives in
`src/python/meshioplusplus/_blender.py`, ships in the wheel, and is tested
without Blender in the default CI matrix. That split is deliberate: the
ParaView plugin under `tools/` is broken precisely because it lives outside the
package — pip never installs it and no test imports it.

## Build

```sh
python -m build --wheel --outdir dist               # the meshio++ wheel
python src/blender/build_extension.py \
    --platform linux-x64 --local-wheel 'dist/*.whl' --outdir dist
```

`--platform` is one of `linux-x64`, `windows-x64`, `macos-x64`, `macos-arm64`.
Without `--local-wheel` the script downloads the matching wheel from PyPI
instead — use the local one on any branch that bumps the version, since the
wheel the manifest must name does not exist on PyPI yet.

Install the resulting zip through `Edit > Preferences > Add-ons > Install from
Disk`, or `blender --command extension install-file dist/*.zip`.

## Develop against a source checkout

Wheel packaging and add-on logic fail in different ways, so separate them:
install meshio++ into Blender's own Python once, and then iterate on `addon/`
with no zip in the loop.

```sh
"$BLENDER/python/bin/python3.11" -m pip install -e /path/to/meshioplusplus
ln -s "$PWD/src/blender/addon" ~/.config/blender/4.2/scripts/addons/meshioplusplus
```

## What is generated, and why

`blender_manifest.toml`, `wheels/` and `dist/` are all gitignored. The manifest
is generated rather than committed because the bundled wheel **file names**
carry the meshio++ version, so a hand-maintained `wheels = [...]` would rot on
the next version bump — and `CLAUDE.md`'s version-bump ritual is already ten
files long, with two of them cited as having silently drifted. The
`blender-packaging` CI job asserts that the generated manifest's version equals
`pyproject.toml`'s and that its wheel list matches the zip exactly, so the
invariant is machine-checked instead of remembered.

## Which wheels are bundled

`meshioplusplus` (cp311, per platform) plus `rich` and its own dependencies.

`rich` is not optional: `meshioplusplus/__init__.py` imports `_cli` eagerly and
`_cli` imports rich. `numpy` is deliberately **not** bundled — Blender ships
it, and a second numpy ahead of Blender's own is a well-known way to break
`bpy`'s numpy interop.

## Tests

The bridge's own tests are `tests/python/test_blender.py`; the pure half runs
everywhere, the `bpy` half behind `pytest.importorskip`. There are no tests
under this directory — `addon/` is option-gathering and menu registration, and
what it calls is tested where it lives.
