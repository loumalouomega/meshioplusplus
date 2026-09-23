"""Freeze pyNastran's reading of the MSC Nastran HDF5 fixtures.

Writes ``tests/python/meshes/nastran_h5/pynastran_reference.npz``: for each nodal
translation/rotation table (EIGENVECTOR, DISPLACEMENT, SPC_FORCE, ...) of each
domain, the node ids and the six components pyNastran's own ``.h5`` reader
(``pyNastran/dev/h5``) returns. ``test_nastran_h5.py`` compares meshio++'s arrays
against it, so the cross-check runs without pyNastran installed.

pyNastran is not a dependency. Run this in a throw-away environment:

    uv venv --python 3.12 /tmp/pyn && uv pip install --python /tmp/pyn/bin/python \
        pyNastran==1.4.1 h5py vtk pandas
    # The wheel has no pyNastran.dev: copy pyNastran/dev/h5 (and dev/__init__.py)
    # from the GitHub repository into the venv's site-packages/pyNastran/dev.
    /tmp/pyn/bin/python tools/gen_nastran_h5_reference.py

Two things in that development reader need working around, and both are done
here rather than by editing it:

- ``_load_results`` has its nodal branch commented out, so
  ``_load_nodal_results`` is called directly.
- ``load_eigenvector`` asserts that every row of a table is a GRID, which the
  SPOINT rows MSC writes into ``EIGENVECTOR`` break. The assertion only guards
  its own geometry wiring (the values are sliced from the raw rows either way),
  so it is disabled for the run by rewriting the one line in memory.
"""

import inspect
import pathlib
import textwrap

import h5py
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
FIXTURES = HERE.parent / "tests" / "python" / "meshes" / "nastran_h5"
FILES = [
    "static_elements",
    "modes_elements",
    "buckling_solid_shell_bar",
    "time_thermal_elements",
]
NAMES = {
    "Eigenvector": "EIGENVECTOR",
    "Displacement": "DISPLACEMENT",
    "SPC Force": "SPC_FORCE",
    "MPC Force": "MPC_FORCE",
    "Applied Load": "APPLIED_LOAD",
    "Velocity": "VELOCITY",
}


def _patched_reader():
    from pyNastran.dev.h5 import read_h5

    src = inspect.getsource(read_h5.load_eigenvector)
    src = src.replace("assert len(_nids) == len(ids)", "pass")
    namespace = vars(read_h5)
    exec(textwrap.dedent(src), namespace)  # noqa: S102 - rebinds load_eigenvector
    return read_h5


def main():
    from pyNastran.dev.h5.h5_utils import h5py_to_dataframe
    from pyNastran.op2.op2 import OP2

    read_h5 = _patched_reader()
    out = {}
    for name in FILES:
        reader = read_h5.pyNastranH5(add_aero=False, add_constraints=False)
        f = h5py.File(FIXTURES / f"{name}.h5", "r")
        reader._load_geometry(f)
        geom = reader.geom_model
        results = {}
        reader._load_nodal_results(
            0,
            results,
            geom,
            OP2(),
            geom._node_ids,
            h5py_to_dataframe(f["/NASTRAN/RESULT/DOMAINS"]),
            f,
        )
        for v in results.values():
            base = v.name.split(" (")[0].split(":")[0]
            if base not in NAMES or not hasattr(v, "TX"):
                continue
            key = f"{name}/{NAMES[base]}/{int(v.domain)}"
            out[key + "/ids"] = np.asarray(v.ID, dtype=np.int64)
            out[key + "/values"] = np.stack(
                [
                    np.asarray(getattr(v, c))
                    for c in ("TX", "TY", "TZ", "RX", "RY", "RZ")
                ],
                axis=1,
            )
    np.savez_compressed(FIXTURES / "pynastran_reference.npz", **out)
    print(f"wrote {len(out) // 2} tables to {FIXTURES / 'pynastran_reference.npz'}")


if __name__ == "__main__":
    main()
