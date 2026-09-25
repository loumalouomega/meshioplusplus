"""Write the LS-DYNA binout fixtures and freeze lasso-python's reading of them.

``tests/python/meshes/lsdyna_binout/`` holds

* ``binout_glstat``: a binout LS-DYNA wrote (``glstat`` and ``rwforc``), from
  Ansys' example data (``result_files/binout_glstat``, MIT), copied as it is;
* ``nodout/binout``: the first 60 outputs of ``nodout`` and ``elout/beam`` and
  the first 7 of ``glstat`` and ``matsum`` of that data's ``binout_matsum``
  (15 MB, a Hybrid III dummy), rewritten with lasso-python's LSDA writer.

``binout_reference.npz`` freezes lasso-python's ``Binout`` reading of both.
Run it in a throwaway environment with ``pip install lasso-python``::

    python tools/gen_binout_reference.py <ansys example-data checkout>/result_files
"""

import os
import shutil
import sys

import numpy as np
from lasso.dyna import Binout
from lasso.dyna.lsda_py3 import Lsda

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "tests", "python", "meshes", "lsdyna_binout")

# (database path, outputs kept)
_SUBSET = (
    (("nodout",), 60),
    (("elout", "beam"), 60),
    (("glstat",), 7),
    (("matsum",), 7),
)


def _children(symbol):
    return {
        (k.decode() if isinstance(k, bytes) else k): v
        for k, v in symbol.children.items()
    }


def _copy_folder(_unused, dst, folder):
    """Every variable of ``folder`` (a Symbol) into ``dst``'s current directory."""
    for name, var in sorted(_children(folder).items()):
        if var.type == 0:
            continue
        dst.write(name, var.type, var.read())


def write_subset(src_path, dst_path):
    src = Lsda(src_path)
    if os.path.exists(dst_path):
        os.remove(dst_path)
    dst = Lsda(dst_path, "w")
    root = _children(src.root)
    for path, count in _SUBSET:
        node = root[path[0]]
        for part in path[1:]:
            node = _children(node)[part]
        kids = _children(node)
        base = "/" + "/".join(path)
        dst.cd(base + "/metadata")
        _copy_folder(None, dst, kids["metadata"])
        steps = sorted(k for k in kids if k[:1] == "d" and k[1:].isdigit())[:count]
        for step in steps:
            dst.cd(f"{base}/{step}")
            _copy_folder(None, dst, kids[step])
    dst.close()


def freeze():
    out = {}
    for key, path in (
        ("glstat", os.path.join(OUT, "binout_glstat")),
        ("nodout", os.path.join(OUT, "nodout", "binout")),
    ):
        b = Binout(path)
        for db in b.read():
            for var in b.read(db):
                try:
                    value = b.read(db, var)
                except Exception:  # noqa: BLE001 - nested databases (elout/beam)
                    continue
                if isinstance(value, list):
                    for sub in value:
                        try:
                            v = np.asarray(b.read(db, var, sub))
                        except Exception:  # noqa: BLE001
                            continue
                        if v.dtype.kind in "iuf":
                            out[f"{key}|{db}/{var}|{sub}"] = v
                    continue
                v = np.asarray(value)
                if v.dtype.kind in "iuf":
                    out[f"{key}|{db}|{var}"] = v
    np.savez_compressed(os.path.join(OUT, "binout_reference.npz"), **out)
    print(len(out), "arrays")


def main():
    if len(sys.argv) > 1:
        src = sys.argv[1]
        shutil.copyfile(
            os.path.join(src, "binout_glstat"), os.path.join(OUT, "binout_glstat")
        )
        os.makedirs(os.path.join(OUT, "nodout"), exist_ok=True)
        write_subset(
            os.path.join(src, "binout_matsum"), os.path.join(OUT, "nodout", "binout")
        )
    freeze()


if __name__ == "__main__":
    main()
