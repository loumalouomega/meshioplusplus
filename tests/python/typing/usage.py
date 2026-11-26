"""Typed-usage probe (roadmap 3.7.4): reads, cleans, and writes a mesh."""

from pathlib import Path

import numpy as np

import meshioplusplus as mp


def roundtrip_clean(src: Path, dst: Path) -> int:
    mesh: mp.Mesh = mp.read(src)
    pts = np.asarray(mesh.points, dtype=np.float64)
    assert pts.ndim == 2
    cleaned: mp.Mesh = mp.clean(mesh)
    mp.write(dst, cleaned)
    return int(pts.shape[0])
