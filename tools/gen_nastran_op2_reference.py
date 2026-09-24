"""Copy the Nastran OP2 fixtures and freeze pyNastran's reading of them.

The OP2 fixtures under ``tests/python/meshes/nastran_op2/`` are pyNastran's own
test models (BSD-3; real MSC and NX output, 32- and 64-bit), copied as they are,
plus ``solid_bending_no_geom.op2``: ``solid_bending.op2`` with its geometry
tables (``GEOM*``, ``EPT*``) removed, next to a copy of its input deck, for the
reader's sibling-deck route.

For every fixture the script reads each result with pyNastran and stores what it
reads in ``pynastran_reference.npz``, under meshio++'s names: per step (subcase,
analysis code, mode or time) the nodal vectors by GRID id and the centre values
of the element stress and strain by element id. The tests compare meshio++
against that file, so they need neither pyNastran nor Nastran. Run it in a
throwaway environment with ``pip install pyNastran==1.4.1``::

    python tools/gen_nastran_op2_reference.py <pyNastran checkout>/models
"""

import os
import shutil
import struct
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "tests", "python", "meshes", "nastran_op2")

FIXTURES = {
    # name in the fixture folder: path under pyNastran's models/
    "static_solid_shell_bar.op2": "sol_101_elements/static_solid_shell_bar.op2",
    "mode_solid_shell_bar.op2": "sol_101_elements/mode_solid_shell_bar.op2",
    "buckling_solid_shell_bar.op2": "sol_101_elements/buckling_solid_shell_bar.op2",
    "solid_bending.op2": "solid_bending/solid_bending.op2",
    "solid_bending.bdf": "solid_bending/solid_bending.bdf",
    "d173.op2": "msc/64_bit/d173.op2",
    "sol401_tstep1.op2": "nx/sol401/sol401_tstep1.op2",
    "time_thermal_elements.op2": "elements/time_thermal_elements.op2",
    "time_thermal_elements_sort2_nx.op2": "elements/time_thermal_elements_sort2_nx.op2",
    "static_elements.op2": "elements/static_elements.op2",
    "ctetra10.op2": "unit/pload4/ctetra10.op2",
}

_NODAL = {
    "displacements": "DISPLACEMENT",
    "eigenvectors": "EIGENVECTOR",
    "spc_forces": "SPC_FORCE",
    "mpc_forces": "MPC_FORCE",
    "load_vectors": "APPLIED_LOAD",
    "temperatures": "TEMPERATURE",
}
_SOLID = {
    "oxx": "X",
    "oyy": "Y",
    "ozz": "Z",
    "txy": "TXY",
    "tyz": "TYZ",
    "txz": "TZX",
    "exx": "X",
    "eyy": "Y",
    "ezz": "Z",
    "exy": "TXY",
    "eyz": "TYZ",
    "exz": "TZX",
    "von_mises": "VON_MISES",
}
_PLATE = {
    "fiber_distance": "FD",
    "fiber_curvature": "FD",
    "oxx": "X",
    "oyy": "Y",
    "txy": "TXY",
    "angle": "ANGLE",
    "omax": "MAJOR",
    "omin": "MINOR",
    "von_mises": "VON_MISES",
    "max_shear": "MAX_SHEAR",
    "exx": "X",
    "eyy": "Y",
    "exy": "TXY",
    "emax": "MAJOR",
    "emin": "MINOR",
    "evm": "VON_MISES",
}
_ROD = {"axial": "A", "SMa": "MSA", "torsion": "T", "SMt": "MST"}
_BAR = dict(
    zip(
        ["s1a", "s2a", "s3a", "s4a", "axial", "smaxa", "smina", "MS_tension"]
        + ["s1b", "s2b", "s3b", "s4b", "smaxb", "sminb", "MS_compression"],
        ["X1A", "X2A", "X3A", "X4A", "AX", "MAXA", "MINA", "MST"]
        + ["X1B", "X2B", "X3B", "X4B", "MAXB", "MINB", "MSC"],
    )
)


def _strip_geometry(src, dst):
    """Copy an OP2 without its GEOM*/EPT* tables (4-byte-word files only)."""
    data = open(src, "rb").read()
    blocks = []
    pos = 0
    while pos < len(data):
        (n,) = struct.unpack("<i", data[pos : pos + 4])
        blocks.append((pos, pos + 8 + n, data[pos + 4 : pos + 4 + n]))
        pos += 8 + n

    def marker(i):
        payload = blocks[i][2]
        return struct.unpack("<i", payload)[0] if len(payload) == 4 else None

    keep = []
    i = 0
    # header: [3] date [7] tape [n] version... [-1] [0]
    while not (marker(i) == -1 and marker(i + 1) == 0):
        keep.append(i)
        i += 1
    keep += [i, i + 1]
    i += 2
    while i < len(blocks):
        if marker(i) == 0:
            keep.append(i)
            i += 1
            continue
        start = i
        name = blocks[i + 1][2][:8].decode("latin-1").strip()
        i += 2
        while (
            (blocks[i][1] - blocks[i][0]) == 12
            and marker(i) is not None
            and marker(i) > 0
        ):
            i += 2  # continuation blocks of the name record
        # the reader's grammar: 0 ends the table, -1 is alone, other negative
        # markers come in triples, a positive one opens a record
        while True:
            m = marker(i)
            if m == 0:
                i += 1
                break
            if m == -1:
                i += 1
            elif m < 0:
                i += 3
            else:
                i += 2
                while marker(i) is not None and marker(i) > 0:
                    i += 2
        if not name.startswith(("GEOM", "EPT")):
            keep.extend(range(start, i))
    with open(dst, "wb") as f:
        for k in keep:
            f.write(data[blocks[k][0] : blocks[k][1]])


def copy_fixtures(models):
    os.makedirs(OUT, exist_ok=True)
    for name, rel in FIXTURES.items():
        shutil.copyfile(os.path.join(models, rel), os.path.join(OUT, name))
    _strip_geometry(
        os.path.join(OUT, "solid_bending.op2"),
        os.path.join(OUT, "solid_bending_no_geom.op2"),
    )
    shutil.copyfile(
        os.path.join(OUT, "solid_bending.bdf"),
        os.path.join(OUT, "solid_bending_no_geom.bdf"),
    )


def _key(key):
    """(subcase, analysis code or None) of a pyNastran result key."""
    if isinstance(key, tuple):
        return int(key[0]), int(key[1])
    return int(key), None


def _tag(obj, ana, it):
    """The step tag: the mode number of a mode (analysis 2 and 9; post-buckling's
    load step, analysis 8), else the index of the time or load step."""
    if ana in (2, 9) and getattr(obj, "modes", None) is not None:
        return int(obj.modes[it])
    if ana == 8 and getattr(obj, "lsdvmns", None) is not None:
        return int(obj.lsdvmns[it])
    return it


def freeze():
    from pyNastran.op2.op2_geom import read_op2_geom

    warnings.filterwarnings("ignore")
    out = {}
    for name in sorted(list(FIXTURES) + ["solid_bending_no_geom.op2"]):
        if not name.endswith(".op2"):
            continue
        path = os.path.join(OUT, name)
        try:
            model = read_op2_geom(path, debug=None, log=None)
        except Exception:
            model = read_op2_geom(path, debug=None, log=None, xref=False)
        stem = name[:-4]
        for attr, label in _NODAL.items():
            for key, obj in getattr(model, attr, {}).items():
                sub, ana = _key(key)
                ana = int(obj.analysis_code) if ana is None else ana
                for it, _t in enumerate(np.atleast_1d(obj._times)):
                    tag = _tag(obj, ana, it)
                    base = f"{stem}|{sub}|{ana}|{tag}|{label}"
                    out[base + "|ids"] = obj.node_gridtype[:, 0]
                    out[base] = (
                        obj.data[it, :, :3]
                        if label != "TEMPERATURE"
                        else obj.data[it, :, 0]
                    )
        for group in ("stress", "strain"):
            container = getattr(model.op2_results, group)
            for attr in sorted(dir(container)):
                if attr.startswith("_") or not attr.endswith("_" + group):
                    continue
                results = getattr(container, attr)
                if not isinstance(results, dict):
                    continue
                etype = attr.split("_")[0]
                for key, obj in results.items():
                    if not hasattr(obj, "data"):
                        continue
                    sub, ana = _key(key)
                    ana = int(obj.analysis_code) if ana is None else ana
                    heads = [str(h) for h in obj.get_headers()]
                    if hasattr(obj, "element_node"):
                        en = obj.element_node
                        rows = np.where(en[:, 1] == 0)[0]
                        eids = en[rows, 0]
                    elif hasattr(obj, "element"):
                        rows = np.arange(len(obj.element))
                        eids = obj.element
                    else:
                        continue
                    if etype in ("ctetra", "chexa", "cpenta", "cpyram"):
                        mapping = _SOLID
                        fibers = [None] * len(rows)
                    elif etype in (
                        "cquad4",
                        "ctria3",
                        "cquad8",
                        "ctria6",
                        "ctriar",
                        "cquadr",
                    ):
                        mapping = _PLATE
                        seen = {}
                        fibers = []
                        for e in eids.tolist():
                            seen[e] = seen.get(e, 0) + 1
                            fibers.append(seen[e])
                    elif etype in ("crod", "conrod"):
                        mapping, fibers = _ROD, [None] * len(rows)
                    elif etype == "ctube":
                        mapping = {
                            "axial": "AS",
                            "SMa": "MSA",
                            "torsion": "TS",
                            "SMt": "MST",
                        }
                        fibers = [None] * len(rows)
                    elif etype == "cbar":
                        mapping, fibers = _BAR, [None] * len(rows)
                    else:
                        continue
                    for it, _t in enumerate(np.atleast_1d(obj._times)):
                        tag = _tag(obj, ana, it)
                        for c, h in enumerate(heads):
                            member = mapping.get(h)
                            if member is None:
                                continue
                            by_name = {}
                            for r, e, fib in zip(rows.tolist(), eids.tolist(), fibers):
                                nm = f"{group.upper()}:{member}{fib if fib else ''}"
                                by_name.setdefault(nm, ([], []))
                                by_name[nm][0].append(e)
                                by_name[nm][1].append(obj.data[it, r, c])
                            for nm, (es, vs) in by_name.items():
                                base = f"{stem}|{sub}|{ana}|{tag}|{nm}|{etype}"
                                out[base + "|ids"] = np.asarray(es, dtype=np.int64)
                                out[base] = np.asarray(vs, dtype=np.float64)
        print(stem, "done")
    np.savez_compressed(os.path.join(OUT, "pynastran_reference.npz"), **out)


def main():
    if len(sys.argv) > 1:
        copy_fixtures(sys.argv[1])
    freeze()


if __name__ == "__main__":
    main()
