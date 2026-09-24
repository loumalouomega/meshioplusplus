"""Freeze pyNastran's coordinate-system transforms for the Nastran result readers.

No MSC or NX fixture puts a GRID in a local system, so this derives two probes
from the committed files and records what pyNastran makes of them:

- ``tests/python/meshes/nastran_h5/cord_reference.npz``: the systems (a tilted
  CORD2C, a CORD2S defined in that cylindrical system, a CORD1R through three
  GRIDs, one of them itself in the CORD2C), the GRID coordinates, CP and CD to
  write into a copy of ``static_elements.h5`` (``test_nastran_h5.py`` does the
  writing), and pyNastran's basic positions and displacements.
- ``tests/python/meshes/nastran_op2/static_solid_shell_bar_cord.op2``: that
  file with every GRID's CP and CD words set to its own CORD2R/C/S systems,
  and ``cord_reference.npz`` beside it with pyNastran's basic positions and
  displacements of subcase 1.

The positions are pyNastran 1.4.1's (``GRID.get_position``). Its displacement
transform composes the cylindrical and spherical rotations the wrong way round
(fixed upstream since); the displacements here use its systems (``beta``,
``origin``) with the rotation matrices of the upstream
``pyNastran/femutils/coord_transforms.py``: ``v = beta^T @ R @ v_local``.

pyNastran is not a dependency. Run this in a throw-away environment with
pyNastran 1.4.1 and h5py, pointing ``PYNASTRAN_SRC`` at a checkout of the
pyNastran repository (for ``femutils/coord_transforms.py``)::

    PYNASTRAN_SRC=.../pyNastran python tools/gen_nastran_cord_reference.py
"""

import importlib.util
import os
import pathlib
import shutil
import struct
import tempfile

import h5py
import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2_geom import read_op2_geom

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE.parent / "tests" / "python" / "meshes"

CORD2C = (5, 0, [1.0, 2.0, 3.0], [1.0, 3.0, 4.0], [2.0, 2.0, 3.0])
CORD2S = (6, 5, [1.0, 30.0, 0.5], [1.0, 30.0, 2.5], [3.0, 60.0, 0.5])


def _upstream_rotations():
    path = (
        pathlib.Path(os.environ["PYNASTRAN_SRC"])
        / "pyNastran"
        / "femutils"
        / "coord_transforms.py"
    )
    spec = importlib.util.spec_from_file_location("coord_transforms", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _to_basic(coords, cd, xyz_basic, values, ct):
    """``v = beta^T @ R @ v_local`` per node, translations and rotations."""
    out = values.copy()
    for i, c in enumerate(cd):
        if c == 0:
            continue
        coord = coords[int(c)]
        beta = coord.beta()
        local = coord.xyz_to_coord_array(
            ((xyz_basic[i] - coord.origin) @ beta.T)[None]
        )[0]
        if coord.type in ("CORD2C", "CORD1C"):
            rot = ct.cylindrical_rotation_matrix(np.radians([local[1]]))[0]
        elif coord.type in ("CORD2S", "CORD1S"):
            rot = ct.spherical_rotation_matrix(
                np.radians([local[1]]), np.radians([local[2]])
            )[0]
        else:
            rot = np.eye(3)
        m = beta.T @ rot
        out[i, :3] = m @ values[i, :3]
        out[i, 3:] = m @ values[i, 3:]
    return out


def hdf5_probe(ct):
    src = MESHES / "nastran_h5" / "static_elements.h5"
    with h5py.File(src, "r") as f:
        grid = f["NASTRAN/INPUT/NODE/GRID"][()]
        disp = f["NASTRAN/RESULT/NODAL/DISPLACEMENT"][()]
    n = len(grid)
    ids = grid["ID"].astype(np.int64)
    x = grid["X"].astype(np.float64).copy()
    cp = np.zeros(n, np.int64)
    cd = np.zeros(n, np.int64)
    for i in range(n):
        k = 1 if i == 1 else i % 4
        if k == 1:
            cp[i] = 5
            x[i] = [1.0 + abs(x[i, 0]), 25.0 * i % 360, x[i, 2]]
        elif k == 2 and i > 2:
            cp[i] = 6
            x[i] = [0.5 + abs(x[i, 1]), 10.0 + 7.0 * i % 170, 15.0 * i % 360]
        elif k == 3 and i > 2:
            cp[i] = 7
        cd[i] = [0, 5, 6, 7][(i + 1) % 4]
    g123 = [int(g) for g in ids[:3]]
    model = BDF(debug=False)
    model.add_cord2c(CORD2C[0], CORD2C[2], CORD2C[3], CORD2C[4], rid=CORD2C[1])
    model.add_cord2s(CORD2S[0], CORD2S[2], CORD2S[3], CORD2S[4], rid=CORD2S[1])
    model.add_cord1r(7, *g123)
    for i in range(n):
        model.add_grid(int(ids[i]), x[i].tolist(), cp=int(cp[i]), cd=int(cd[i]))
    model.cross_reference()
    xyz0 = np.array([model.nodes[int(g)].get_position() for g in ids])
    first = disp["DOMAIN_ID"][0]
    pos = {int(g): i for i, g in enumerate(ids)}
    local = np.full((n, 6), np.nan)
    for r in disp[disp["DOMAIN_ID"] == first]:
        if int(r["ID"]) in pos:
            local[pos[int(r["ID"])]] = [r[m] for m in ("X", "Y", "Z", "RX", "RY", "RZ")]
    basic = _to_basic(model.coords, cd, xyz0, local, ct)
    np.savez(
        MESHES / "nastran_h5" / "cord_reference.npz",
        x_local=x,
        cp=cp,
        cd=cd,
        cord2=np.array([[c[0], c[1], *c[2], *c[3], *c[4]] for c in (CORD2C, CORD2S)]),
        cord2_type=np.array([2, 3]),
        cord1r=np.array([[7, *g123]]),
        xyz_basic=xyz0,
        displacement_basic=basic,
    )


def op2_probe(ct):
    src = MESHES / "nastran_op2" / "static_solid_shell_bar.op2"
    out = MESHES / "nastran_op2" / "static_solid_shell_bar_cord.op2"
    data = bytearray(src.read_bytes())
    # The GRID record: key (4501, 45, 1) behind its Fortran length word.
    key = struct.pack("<3i", 4501, 45, 1)
    at = bytes(data).find(key)
    assert at > 0 and bytes(data).find(key, at + 1) < 0
    size = struct.unpack("<i", data[at - 4 : at])[0]
    words = np.frombuffer(bytes(data[at : at + size]), "<i4").copy()
    body = words[3:].reshape(-1, 8)
    cids = [0, 12, 13, 11, 2, 3]  # the file's CORD2C/2S/2R systems
    for i in range(len(body)):
        body[i, 1] = cids[i % 6]
        body[i, 5] = cids[(i + 2) % 6]
    data[at : at + size] = words.tobytes()
    out.write_bytes(bytes(data))
    with tempfile.TemporaryDirectory() as tmp:
        copy = pathlib.Path(tmp) / out.name
        shutil.copyfile(out, copy)
        model = read_op2_geom(str(copy), debug=None, log=None)
    nids = np.array(sorted(model.nodes), dtype=np.int64)
    xyz0 = np.array([model.nodes[int(g)].get_position() for g in nids])
    result = model.displacements[sorted(model.displacements)[0]]
    dn = result.node_gridtype[:, 0].astype(np.int64)
    values = result.data[0].astype(np.float64)
    keep = np.isin(dn, nids)
    rows = np.searchsorted(nids, dn[keep])
    cd = np.array([model.nodes[int(g)].cd for g in dn[keep]])
    basic = _to_basic(model.coords, cd, xyz0[rows], values[keep], ct)
    np.savez(
        MESHES / "nastran_op2" / "cord_reference.npz",
        nids=nids,
        xyz_basic=xyz0,
        displacement_nids=dn[keep],
        displacement_basic=basic,
    )


def main():
    ct = _upstream_rotations()
    hdf5_probe(ct)
    op2_probe(ct)


if __name__ == "__main__":
    main()
