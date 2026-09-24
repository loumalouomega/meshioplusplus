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
    "freq_elements2.op2": "elements/freq_elements2.op2",
    "modes_complex_elements.op2": "elements/modes_complex_elements.op2",
    "test_vba.op2": "nx/test_vba/test_vba.op2",
}
# Random result groups of pyNastran's op2_results -> meshio++'s suffix.
_RANDOM = {"psd": "_PSD", "ato": "_ATO", "rms": "_RMS", "no": "_NO", "crm": "_CRM"}

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
    "ovm": "VON_MISES",
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


_COMPOSITE = dict(
    zip(
        ["o11", "o22", "t12", "t1z", "t2z", "angle", "major", "minor", "max_shear"]
        + ["e11", "e22", "e12", "e1z", "e2z", "von_mises"],
        ["X1", "Y1", "T1", "L1", "L2", "ANGLE", "MAJOR", "MINOR", "MAX_SHEAR"]
        + ["X1", "Y1", "T1", "L1", "L2", "VON_MISES"],
    )
)
_BEAM = dict(
    zip(
        ["sxc", "sxd", "sxe", "sxf", "smax", "smin", "MS_tension", "MS_compression"],
        ["XC", "XD", "XE", "XF", "MAX", "MIN", "MST", "MSC"],
    )
)
_GPF = dict(
    zip(["f1", "f2", "f3", "m1", "m2", "m3"], ["F1", "F2", "F3", "M1", "M2", "M3"])
)


def _put(out, base, ids, cols, values):
    out[base + "|ids"] = np.asarray(ids, dtype=np.int64)
    out[base + "|cols"] = np.asarray(cols, dtype=np.float64)
    out[base] = np.asarray(values, dtype=np.float64)


def _wide(model, stem, out):
    """The multi-valued results: composite plies (``|cols`` = ply id), plate and
    solid corners (``|cols`` = GRID id), CBEAM stations (``|cols`` = station
    distance) and grid point forces (``|cols`` = GRID id; element 0 rows by GRID
    under ``GRID_FORCE:<label>``)."""
    for group in ("stress", "strain"):
        container = getattr(model.op2_results, group)
        for attr in sorted(dir(container)):
            if attr.startswith("_") or not attr.endswith("_" + group):
                continue
            results = getattr(container, attr)
            if not isinstance(results, dict):
                continue
            etype = attr.split("_")[0]
            composite = "_composite_" in attr
            for key, obj in results.items():
                if not hasattr(obj, "data") or getattr(obj, "is_complex", False):
                    continue
                sub, ana = _key(key)
                ana = int(obj.analysis_code) if ana is None else ana
                heads = [str(h) for h in obj.get_headers()]
                if composite and hasattr(obj, "element_layer"):
                    mapping, suffix = _COMPOSITE, "@ply"
                    eids, cols = obj.element_layer[:, 0], obj.element_layer[:, 1]
                    fibers = [None] * len(eids)
                elif etype == "cbeam" and hasattr(obj, "xxb"):
                    mapping, suffix = _BEAM, "@station"
                    eids, cols = obj.element_node[:, 0], obj.xxb
                    fibers = [None] * len(eids)
                elif hasattr(obj, "element_node") and etype in (
                    "ctetra",
                    "chexa",
                    "cpenta",
                    "cpyram",
                    "cquad4",
                    "cquad8",
                    "ctria6",
                    "ctriar",
                    "cquadr",
                ):
                    en = obj.element_node
                    plate = etype not in ("ctetra", "chexa", "cpenta", "cpyram")
                    mapping, suffix = (_PLATE if plate else _SOLID), "@corner"
                    keep = np.where(en[:, 1] != 0)[0]
                    if len(keep) == 0:
                        continue
                    eids, cols = en[keep, 0], en[keep, 1]
                    if plate:
                        seen = {}
                        fibers = []
                        for e, n in zip(eids.tolist(), cols.tolist()):
                            seen[(e, n)] = seen.get((e, n), 0) + 1
                            fibers.append(seen[(e, n)])
                    else:
                        fibers = [None] * len(eids)
                    rows = keep
                else:
                    continue
                if suffix != "@corner":
                    rows = np.arange(len(eids))
                for it, _t in enumerate(np.atleast_1d(obj._times)):
                    tag = _tag(obj, ana, it)
                    for c, h in enumerate(heads):
                        member = mapping.get(h)
                        if member is None:
                            continue
                        by_name = {}
                        for r, e, col, fib in zip(
                            rows.tolist(), eids.tolist(), cols.tolist(), fibers
                        ):
                            nm = f"{group.upper()}:{member}{fib if fib else ''}{suffix}"
                            entry = by_name.setdefault(nm, ([], [], []))
                            entry[0].append(e)
                            entry[1].append(col)
                            entry[2].append(obj.data[it, r, c])
                        for nm, (es, cs, vs) in by_name.items():
                            _put(
                                out,
                                f"{stem}|{sub}|{ana}|{tag}|{nm}|{etype}",
                                es,
                                cs,
                                vs,
                            )
                    if etype == "cbeam" and suffix == "@station":
                        _put(
                            out,
                            f"{stem}|{sub}|{ana}|{_tag(obj, ana, it)}|{group.upper()}:SD@station|{etype}",
                            eids,
                            cols,
                            cols,
                        )
    for key, obj in getattr(model, "grid_point_forces", {}).items():
        if np.iscomplexobj(obj.data):  # complex grid point forces are not read
            continue
        sub, ana = _key(key)
        ana = int(obj.analysis_code) if ana is None else ana
        ne = obj.node_element
        names = [str(n).replace(" ", "").replace("*", "") for n in obj.element_names]
        for it, _t in enumerate(np.atleast_1d(obj._times)):
            tag = _tag(obj, ana, it)
            ne_it = ne[it] if ne.ndim == 3 else ne
            names_it = (
                names
                if np.ndim(obj.element_names) == 1
                else [
                    str(n).replace(" ", "").replace("*", "")
                    for n in obj.element_names[it]
                ]
            )
            for c, h in enumerate(obj.get_headers()):
                member = _GPF.get(str(h))
                if member is None:
                    continue
                elem = {}
                other = {}
                for r, (nid, eid) in enumerate(ne_it.tolist()):
                    v = obj.data[it, r, c]
                    if eid > 0:
                        entry = elem.setdefault(f"GRID_FORCE:{member}", ([], [], []))
                        entry[0].append(eid)
                        entry[1].append(nid)
                        entry[2].append(v)
                    else:
                        entry = other.setdefault(
                            f"GRID_FORCE:{names_it[r]}:{member}", ([], [], [])
                        )
                        entry[0].append(nid)
                        entry[1].append(nid)
                        entry[2].append(v)
                for nm, (es, cs, vs) in list(elem.items()) + list(other.items()):
                    _put(out, f"{stem}|{sub}|{ana}|{tag}|{nm}|gpf", es, cs, vs)


# Element forces (``force.<type>_force``), heat fluxes
# (``thermal_load.<type>_thermal_load``) and energies (``strain_energy.<type>_
# strain_energy``): pyNastran's headers -> meshio++'s members, per group.
_FORCE = {
    "axial": "AF",
    "axial_force": "AF",
    "torsion": "TRQ",
    "torque": "TRQ",
    "spring_force": "F",
    "damper_force": "F",
    "bending_moment_a1": "BM1A",
    "bending_moment_a2": "BM2A",
    "bending_moment_b1": "BM1B",
    "bending_moment_b2": "BM2B",
    "bending_moment_1a": "BM1A",
    "bending_moment_2a": "BM2A",
    "bending_moment_1b": "BM1B",
    "bending_moment_2b": "BM2B",
    "shear1": "TS1",
    "shear2": "TS2",
    "mx": "MX",
    "my": "MY",
    "mxy": "MXY",
    "bmx": "BMX",
    "bmy": "BMY",
    "bmxy": "BMXY",
    "tx": "TX",
    "ty": "TY",
    "fx": "FX",
    "fy": "FY",
    "fz": "FZ",
    "mz": "MZ",
    "sfy": "SFY",
    "sfz": "SFZ",
    "u": "U",
    "v": "V",
    "w": "W",
    "sv": "SV",
    "sw": "SW",
    "force41": "F41",
    "force21": "F21",
    "force12": "F12",
    "force32": "F32",
    "force23": "F23",
    "force43": "F43",
    "force34": "F34",
    "force14": "F14",
    "kick_force1": "KF1",
    "shear12": "S12",
    "kick_force2": "KF2",
    "shear23": "S23",
    "kick_force3": "KF3",
    "shear34": "S34",
    "kick_force4": "KF4",
    "shear41": "S41",
}
_FLUX = dict(
    zip(
        ["xgrad", "ygrad", "zgrad", "xflux", "yflux", "zflux"]
        + ["grad1", "grad2", "grad3", "flux1", "flux2", "flux3"]
        + ["fapplied", "free_conv", "force_conv", "frad", "ftotal"],
        ["XGRAD", "YGRAD", "ZGRAD", "XFLUX", "YFLUX", "ZFLUX"] * 2
        + ["FAPPLIED", "FREECONV", "FORCECONV", "FRAD", "FTOTAL"],
    )
)
_ENERGY = {"strain_energy": "ENERGY", "percent": "PCT", "strain_energy_density": "DEN"}


def _forces(model, stem, out):
    """The centre (or end A, or single) values of the element forces, heat
    fluxes and strain energies by element id, complex ones as ``_real`` and
    ``_imag``. Left out where pyNastran misreads (op2_layouts in the tools'
    notes): complex CSHEAR forces (it pairs the wrong words), CBAR-100 forces
    (they share ``cbar_force`` with CBAR-34)."""
    for prefix, group, mapping in (
        ("force", "ELEMENT_FORCE", _FORCE),
        ("thermal_load", "HEAT_FLUX", _FLUX),
        ("strain_energy", "ENERGY", _ENERGY),
    ):
        container = getattr(model.op2_results, prefix, None)
        if container is None:
            continue
        for attr in sorted(dir(container)):
            if attr.startswith("_"):
                continue
            results = getattr(container, attr)
            if not isinstance(results, dict):
                continue
            etype = attr.split("_")[0]
            for key, obj in results.items():
                if not hasattr(obj, "data") or not hasattr(obj, "get_headers"):
                    continue
                cls = type(obj).__name__
                if "CBar100" in cls or ("Shear" in cls and "Complex" in cls):
                    continue
                sub, ana = _key(key)
                ana = int(obj.analysis_code) if ana is None else ana
                heads = [str(h) for h in obj.get_headers()]
                for it, _t in enumerate(np.atleast_1d(obj._times)):
                    if hasattr(obj, "element_node"):
                        en = obj.element_node
                        first = {}
                        for r, e in enumerate(en[:, 0].tolist()):
                            first.setdefault(e, r)  # the centre, or end A
                        rows = np.array(list(first.values()), dtype=np.int64)
                        eids = np.array(list(first.keys()), dtype=np.int64)
                    elif hasattr(obj, "element"):
                        el = obj.element
                        eids = el[it] if el.ndim == 2 else el
                        rows = np.arange(len(eids))
                        # an energy set's totals row is no element
                        keep = eids != 100000000
                        eids, rows = eids[keep], rows[keep]
                    else:
                        continue
                    tag = _tag(obj, ana, it)
                    for c, h in enumerate(heads):
                        member = mapping.get(h)
                        if member is None or c >= obj.data.shape[2]:
                            continue
                        v = obj.data[it, rows, c]
                        parts = (
                            (("_real", v.real), ("_imag", v.imag))
                            if np.iscomplexobj(v)
                            else (("", v),)
                        )
                        for part, pv in parts:
                            base = f"{stem}|{sub}|{ana}|{tag}|{group}:{member}{part}|{etype}"
                            out[base + "|ids"] = np.asarray(eids, dtype=np.int64)
                            out[base] = np.asarray(pv, dtype=np.float64)


def _strip_geometry(src, dst, prefixes=("GEOM", "EPT")):
    """Copy an OP2 without the tables whose names start with ``prefixes``
    (4-byte-word files only)."""
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
        if not name.startswith(tuple(prefixes)):
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
    # Without GEOM1 (and no deck beside them) the points come from the basic
    # grid point table: BGPDTS + EQEXINS (MSC), and NX's newer BGPDT.
    for name in ("static_elements", "sol401_tstep1"):
        _strip_geometry(
            os.path.join(OUT, name + ".op2"),
            os.path.join(OUT, name + "_bgpdt.op2"),
            ("GEOM1",),
        )
    # The SORT2 tables alone: without their SORT1 twins.
    _strip_geometry(
        os.path.join(OUT, "time_thermal_elements_sort2_nx.op2"),
        os.path.join(OUT, "time_thermal_elements_sort2_only.op2"),
        ("OUGV1", "OPG1", "OEF1X"),
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
    if ana == 8:
        # pyNastran's force objects carry no load steps: meshio++ numbers a
        # post-buckling step from 1, as the steps' own load step words do.
        steps = getattr(obj, "lsdvmns", None)
        return int(steps[it]) if steps is not None else it + 1
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
                    values = (
                        obj.data[it, :, :3]
                        if label != "TEMPERATURE"
                        else obj.data[it, :, 0]
                    )
                    if np.iscomplexobj(values):  # meshio++'s _real and _imag
                        for part, f in (("_real", np.real), ("_imag", np.imag)):
                            out[base + part + "|ids"] = obj.node_gridtype[:, 0]
                            out[base + part] = f(values)
                        continue
                    out[base + "|ids"] = obj.node_gridtype[:, 0]
                    out[base] = values
        # Random results, keyed by their set's time (PSD per frequency; RMS,
        # NO ... once): T=<time> in place of the step tag.
        for kind, suffix in _RANDOM.items():
            group = getattr(model.op2_results, kind, None)
            for attr, label in _NODAL.items():
                results = getattr(group, attr, None) if group is not None else None
                for key, obj in (results or {}).items():
                    sub, ana = _key(key)
                    ana = int(obj.analysis_code) if ana is None else ana
                    for it, t in enumerate(np.atleast_1d(obj._times)):
                        base = f"{stem}|{sub}|{ana}|T={float(t)!r}|{label}{suffix}"
                        out[base + "|ids"] = obj.node_gridtype[:, 0]
                        out[base] = obj.data[it, :, :3]
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
                                v = obj.data[it, r, c]
                                parts = (
                                    (("_real", v.real), ("_imag", v.imag))
                                    if np.iscomplexobj(obj.data)
                                    else (("", v),)
                                )
                                for part, pv in parts:
                                    by_name.setdefault(nm + part, ([], []))
                                    by_name[nm + part][0].append(e)
                                    by_name[nm + part][1].append(pv)
                            for nm, (es, vs) in by_name.items():
                                base = f"{stem}|{sub}|{ana}|{tag}|{nm}|{etype}"
                                out[base + "|ids"] = np.asarray(es, dtype=np.int64)
                                out[base] = np.asarray(vs, dtype=np.float64)
        _forces(model, stem, out)
        _wide(model, stem, out)
        print(stem, "done")
    np.savez_compressed(os.path.join(OUT, "pynastran_reference.npz"), **out)


def main():
    if len(sys.argv) > 1:
        copy_fixtures(sys.argv[1])
    freeze()


if __name__ == "__main__":
    main()
