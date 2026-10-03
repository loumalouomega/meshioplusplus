"""Fresh-process native snapshots: FileSource caches its mmap threshold once."""

import os
import pathlib
import subprocess
import sys

import numpy as np

from meshioplusplus import _core

ROOT = pathlib.Path(__file__).resolve().parents[2]


def mesh_arrays(mesh):
    arrays = {
        "points": mesh.points,
        "cell_types": np.asarray([block.type for block in mesh.cells]),
    }
    for i, block in enumerate(mesh.cells):
        arrays[f"cells/{i}"] = block.data
    for scope in ("point_data", "cell_data", "field_data"):
        for name, data in getattr(mesh, scope).items():
            values = data if scope == "cell_data" else [data]
            for i, value in enumerate(values):
                arrays[f"{scope}/{name}/{i}"] = value
    arrays["regions/meta"] = np.asarray(
        [
            [region.kind, region.name, str(region.dim), str(region.tag)]
            for region in mesh.regions
        ]
    )
    for i, region in enumerate(mesh.regions):
        arrays[f"regions/{i}"] = region.entries
    return arrays


def native_snapshot(path, reader, threshold, tmp_path, extra_sources=()):
    """Read with a threshold set before native initialization, then destroy input."""
    env = os.environ.copy()
    env["MESHIOPLUSPLUS_MMAP_THRESHOLD"] = threshold
    skips = env.get("SKBUILD_EDITABLE_SKIP", "").split(os.pathsep)
    skips.extend(
        str(directory) for directory in (ROOT / "build").glob("*") if directory.is_dir()
    )
    env["SKBUILD_EDITABLE_SKIP"] = os.pathsep.join(skips)
    output = tmp_path / "native-snapshot.npz"
    code = """
import importlib.util, pathlib, sys
spec = importlib.util.spec_from_file_location('meshioplusplus._core', sys.argv[1])
core = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = core
spec.loader.exec_module(core)
from tests.python.text_io_helpers import mesh_arrays
import numpy as np
path = pathlib.Path(sys.argv[2])
result = getattr(core, sys.argv[3])(str(path))
inputs = [path]
if sys.argv[3] == 'ensight_read':
    inputs.append(path.with_suffix('.geo'))
inputs.extend(pathlib.Path(file) for file in sys.argv[5:])
for file in inputs:
    size = file.stat().st_size
    file.unlink()
    file.write_bytes(b'?' * size)
np.savez(sys.argv[4], **mesh_arrays(result))
"""
    subprocess.run(
        [sys.executable, "-c", code, _core.__file__, str(path), reader, str(output)]
        + [str(file) for file in extra_sources],
        cwd=ROOT,
        env=env,
        check=True,
    )
    with np.load(output, allow_pickle=False) as data:
        return {name: data[name] for name in data.files}
