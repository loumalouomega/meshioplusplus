#!/usr/bin/env python3
"""Regenerate the I-DEAS universal-file fixtures under ``tests/python/meshes/unv/``.

Three independent sources, so that the reader is never only checked against its own
writer:

* default: small files written from the dataset record layouts below (SDRL dataset
  pages, as quoted in pyuff's ``get_structure_*`` docstrings);
* ``--gmsh``: quadratic meshes written by gmsh (``pip install gmsh``; 4.15.2 wrote the
  committed files), each as ``<name>.unv`` and its ``<name>.msh`` twin. The node order
  and the groups of the ``.unv`` are checked against the ``.msh``, whose reader is
  already pinned;
* ``--pyuff``: result files written by pyuff (MIT, ``pip install pyuff``; 2.5.6 wrote
  the committed files): units, coordinate systems, legacy nodes, 55 mode shapes and
  58 / 58b frequency-response functions.

    python tools/gen_unv_fixtures.py [--gmsh] [--pyuff]
"""

import argparse
import pathlib

import numpy as np

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "unv"
)


def _nodes(points):
    out = ["    -1", "  2411"]
    for k, p in enumerate(points, 1):
        out.append(f"{k:10d}{1:10d}{1:10d}{11:10d}")
        out.append("".join(f"{x:25.16E}".replace("E", "D") for x in p))
    out.append("    -1")
    return out


def _elements(elements):
    out = ["    -1", "  2412"]
    for label, fe, pid, mid, nodes in elements:
        out.append(f"{label:10d}{fe:10d}{pid:10d}{mid:10d}{7:10d}{len(nodes):10d}")
        if fe in (11, 21, 22, 23, 24, 25):
            out.append(f"{0:10d}{1:10d}{1:10d}")
        for i in range(0, len(nodes), 8):
            out.append("".join(f"{n:10d}" for n in nodes[i : i + 8]))
    out.append("    -1")
    return out


def _group_2467(groups):
    out = ["    -1", "  2467"]
    for number, name, ents in groups:
        out.append(f"{number:10d}" + f"{0:10d}" * 6 + f"{len(ents):10d}")
        out.append(name)
        for i in range(0, len(ents), 2):
            out.append(
                "".join(f"{t:10d}{tag:10d}{0:10d}{0:10d}" for t, tag in ents[i : i + 2])
            )
    out.append("    -1")
    return out


def _group_2429(groups):
    out = ["    -1", "  2429"]
    for number, name, ents in groups:
        out.append(f"{number:10d}" + f"{0:10d}" * 6 + f"{len(ents):10d}")
        out.append(name)
        for i in range(0, len(ents), 4):
            out.append("".join(f"{t:10d}{tag:10d}" for t, tag in ents[i : i + 4]))
    out.append("    -1")
    return out


def _e13(values):
    return "".join(f"{v:13.5E}" for v in values)


def _2414(label, name, location, analysis, char, rtype, ndv, r10, r12, rows):
    out = ["    -1", "  2414", f"{label:10d}", name, f"{location:10d}"]
    out += ["NONE"] * 5
    out.append(f"{1:10d}{analysis:10d}{char:10d}{rtype:10d}{2:10d}{ndv:10d}")
    out.append("".join(f"{v:10d}" for v in r10))
    out.append(f"{0:10d}{0:10d}")
    out.append(_e13(r12[:6]))
    out.append(_e13(r12[6:]))
    for ent, vals in rows:
        out.append(f"{ent:10d}" if location == 1 else f"{ent:10d}{len(vals):10d}")
        for i in range(0, len(vals), 6):
            out.append(_e13(vals[i : i + 6]))
    out.append("    -1")
    return out


# Two unit hexahedra side by side along x.
HEX2_POINTS = [
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    [2, 0, 0], [2, 1, 0], [2, 0, 1], [2, 1, 1],
]  # fmt: skip
HEX2_ELEMENTS = [
    (1, 115, 3, 7, [1, 2, 3, 4, 5, 6, 7, 8]),
    (2, 115, 4, 7, [2, 9, 10, 3, 6, 11, 12, 7]),
    (3, 94, 5, 8, [9, 10, 12, 11]),
]


def modes_2414():
    """Three normal modes of displacement, a static stress on the elements, groups
    in both layouts (2467 quadruples, 2429 pairs) and a mixed node+element group."""
    lines = _nodes(HEX2_POINTS) + _elements(HEX2_ELEMENTS)
    lines += _group_2467(
        [
            (10, "left", [(8, 1)]),
            (11, "right_face", [(7, 9), (7, 10), (7, 11), (7, 12), (8, 3)]),
        ]
    )
    lines += _group_2429([(12, "origin", [(7, 1)])])
    for mode, freq in ((1, 12.5), (2, 31.0), (3, 47.25)):
        r10 = [0, 0, 0, 0, 1, mode, 0, 0]
        r12 = [0.0, freq, (2 * np.pi * freq) ** 2, 1.0, 0.0, 0.0] + [0.0] * 6
        rows = [(n, [0.0, 0.0, mode * 0.1 * n]) for n in range(1, 13)]
        lines += _2414(mode, "Mode shape", 1, 2, 2, 8, 3, r10, r12, rows)
    # Static stress on the two hexahedra (file order Sxx Sxy Syy Sxz Syz Szz).
    rows = [
        (1, [1.0, 4.0, 2.0, 6.0, 5.0, 3.0]),
        (2, [10.0, 40.0, 20.0, 60.0, 50.0, 30.0]),
    ]
    lines += _2414(
        10, "NONE", 2, 1, 4, 2, 6, [0, 0, 0, 0, 1, 0, 0, 0], [0.0] * 12, rows
    )
    return lines


def transient_2414():
    """A temperature at three times, plus two nodes absent from the last step."""
    lines = _nodes(HEX2_POINTS[:8]) + _elements(HEX2_ELEMENTS[:1])
    for step, t in ((1, 0.0), (2, 0.5), (3, 1.0)):
        r10 = [0, 0, 0, 0, 0, 0, step, 0]
        r12 = [t] + [0.0] * 11
        nodes = range(1, 9) if step < 3 else range(1, 7)
        rows = [(n, [100.0 * t + n]) for n in nodes]
        lines += _2414(step, "Temperature", 1, 4, 1, 5, 1, r10, r12, rows)
    return lines


def legacy55():
    """Code-Aster mode of meshio++ <= 15.5: 55/57 with four zero header lines, and
    the element data in 57."""
    lines = _nodes(HEX2_POINTS[:8]) + _elements(HEX2_ELEMENTS[:1])
    lines += ["    -1", "    55"] + ["temp"] * 5
    lines += [
        f"{1:10d}{0:10d}{1:10d}{0:10d}{4:10d}{1:10d}",
        f"{0:10d}" * 8,
        f"{0:10d}" * 2,
    ]
    lines += [_e13([0.0] * 6), _e13([0.0] * 6)]
    for n in range(1, 9):
        lines += [f"{n:10d}", _e13([float(n)])]
    lines += ["    -1", "    -1", "    57"] + ["pressure"] * 5
    lines += [
        f"{1:10d}{0:10d}{1:10d}{0:10d}{4:10d}{1:10d}",
        f"{0:10d}" * 8,
        f"{0:10d}" * 2,
    ]
    lines += [_e13([0.0] * 6), _e13([0.0] * 6), f"{1:10d}", _e13([7.5]), "    -1"]
    return lines


def units_cs_legacy():
    """164 units, a rotated-and-shifted Cartesian 2420 system, 15/780 legacy datasets."""
    lines = ["    -1", "   164", f"{5:10d}{'mm (milli newton)':>20s}{2:10d}"]
    lines += [
        "  1.0000000000000000D+03  1.0000000000000000D+03  1.0000000000000000D+00"
    ]
    lines += ["  2.7315000000000000D+02", "    -1"]
    lines += [
        "    -1",
        "  2420",
        f"{1:10d}",
        "part",
        f"{1:10d}{0:10d}{8:10d}",
        "global",
    ]
    lines += [
        "".join(f"{v:25.16E}" for v in row)
        for row in np.vstack([np.eye(3), np.zeros(3)])
    ]
    lines += [f"{2:10d}{0:10d}{8:10d}", "rotated"]
    rot = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [10.0, 0.0, 0.0]]
    )
    lines += ["".join(f"{v:25.16E}" for v in row) for row in rot]
    lines += ["    -1"]
    # Legacy dataset 15: one line per node, the last node defined in system 2.
    lines += ["    -1", "    15"]
    pts = [(1, 1, [0, 0, 0]), (2, 1, [1, 0, 0]), (3, 1, [0, 1, 0]), (4, 2, [1, 0, 0])]
    for label, cs, x in pts:
        lines.append(f"{label:10d}{cs:10d}{cs:10d}{11:10d}" + _e13(x))
    lines += ["    -1", "    -1", "   780"]
    lines.append("".join(f"{v:10d}" for v in (1, 91, 1, 6, 1, 9, 7, 3)))
    lines.append("".join(f"{v:10d}" for v in (1, 2, 3)))
    lines.append("".join(f"{v:10d}" for v in (2, 21, 1, 6, 1, 9, 7, 2)))
    lines.append("".join(f"{v:10d}" for v in (0, 1, 1)))
    lines.append("".join(f"{v:10d}" for v in (1, 4)))
    lines += ["    -1"]
    return lines


def write_default():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, fn in (
        ("modes_2414.unv", modes_2414),
        ("transient_2414.unv", transient_2414),
        ("legacy55.unv", legacy55),
        ("units_cs_legacy.unv", units_cs_legacy),
    ):
        (OUT / name).write_text("\n".join(fn()) + "\n", newline="\n")
        print("wrote", OUT / name)


def write_gmsh():
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    def one(name, build):
        gmsh.clear()
        gmsh.model.add(name)
        build()
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
        gmsh.option.setNumber("Mesh.SaveGroupsOfNodes", 1)
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(2)
        gmsh.write(str(OUT / f"{name}.unv"))
        gmsh.write(str(OUT / f"{name}.msh"))
        print("wrote", OUT / f"{name}.unv", "and .msh")

    def hexa():
        gmsh.model.occ.addBox(0, 0, 0, 1, 2, 3)
        gmsh.model.occ.synchronize()
        for c in gmsh.model.getEntities(1):
            gmsh.model.mesh.setTransfiniteCurve(c[1], 3)
        for s in gmsh.model.getEntities(2):
            gmsh.model.mesh.setTransfiniteSurface(s[1])
            gmsh.model.mesh.setRecombine(2, s[1])
        gmsh.model.mesh.setTransfiniteVolume(1)
        gmsh.model.addPhysicalGroup(3, [1], 10, "solid")
        gmsh.model.addPhysicalGroup(2, [1], 11, "face")

    def wedge():
        geo = gmsh.model.geo
        for i, (x, y) in enumerate(((0, 0), (1, 0), (0, 1)), 1):
            geo.addPoint(x, y, 0, 1, i)
        for i, (a, b) in enumerate(((1, 2), (2, 3), (3, 1)), 1):
            geo.addLine(a, b, i)
            geo.mesh.setTransfiniteCurve(i, 3)
        geo.addCurveLoop([1, 2, 3], 1)
        geo.addPlaneSurface([1], 1)
        geo.mesh.setTransfiniteSurface(1)
        ex = geo.extrude([(2, 1)], 0, 0, 2, [2], recombine=True)
        geo.synchronize()
        gmsh.model.addPhysicalGroup(3, [ex[1][1]], 20, "prisms")

    def tet():
        gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 0.7)
        gmsh.model.addPhysicalGroup(3, [1], 30, "volume")
        gmsh.model.addPhysicalGroup(2, [1, 2], 31, "sides")

    OUT.mkdir(parents=True, exist_ok=True)
    one("gmsh_hex20", hexa)
    one("gmsh_wedge15", wedge)
    one("gmsh_tet10", tet)
    gmsh.finalize()


def write_pyuff():
    import pyuff

    def fresh(name):
        path = OUT / name
        if path.exists():
            path.unlink()
        return pyuff.UFF(str(path), "w"), path

    nodes = dict(
        node_nums=[1, 2, 3], def_cs=[1, 1, 1], disp_cs=[1, 1, 1], color=[11] * 3
    )
    xyz = dict(x=[0.0, 1.0, 1.0], y=[0.0, 0.0, 1.0], z=[0.0, 0.0, 0.0])

    uff, path = fresh("pyuff_modes.uff")
    uff.write_sets(
        pyuff.prepare_164(
            units_code=1,
            units_description="SI",
            temp_mode=2,
            length=1.0,
            force=1.0,
            temp=1.0,
            temp_offset=273.15,
        ),
        mode="add",
    )
    uff.write_sets(pyuff.prepare_2411(**nodes, **xyz), mode="add")
    for mode, freq in ((1, 10.5), (2, 22.0)):
        uff.write_sets(
            pyuff.prepare_55(
                model_type=1,
                analysis_type=2,
                data_ch=2,
                spec_data_type=8,
                data_type=2,
                n_data_per_node=3,
                node_nums=np.array([1, 2, 3]),
                r1=np.array([1.0, 2.0, 3.0]) * mode,
                r2=np.array([4.0, 5.0, 6.0]) * mode,
                r3=np.array([7.0, 8.0, 9.0]) * mode,
                load_case=1,
                mode_n=mode,
                freq=freq,
                modal_m=1.0,
                modal_damp_vis=0.01,
                modal_damp_his=0.0,
            ),
            mode="add",
        )
    print("wrote", path)

    freq = np.linspace(0.0, 20.0, 9)

    def frfs(binary):
        out = []
        for node in (1, 2, 3):
            for d in (1, 2, 3):
                data = (node + 0.1 * d) * (freq + 1j * (freq + 1))
                out.append(
                    pyuff.prepare_58(
                        binary=binary,
                        id1=f"FRF {node}:{d}",
                        func_type=4,
                        rsp_node=node,
                        rsp_dir=d,
                        ref_node=1,
                        ref_dir=3,
                        abscissa_spacing=1,
                        abscissa_min=0.0,
                        abscissa_inc=2.5,
                        abscissa_spec_data_type=18,
                        ordinate_spec_data_type=12,
                        orddenom_spec_data_type=13,
                        num_pts=len(freq),
                        ord_data_type=6,
                        x=freq,
                        data=data,
                    )
                )
        return out

    for name, binary in (("pyuff_frf.uff", 0), ("pyuff_frf_58b.uff", 1)):
        uff, path = fresh(name)
        uff.write_sets(pyuff.prepare_2411(**nodes, **xyz), mode="add")
        for dset in frfs(binary):
            uff.write_sets(dset, mode="add")
        # A real single-precision coherence on the same grid.
        uff.write_sets(
            pyuff.prepare_58(
                binary=binary,
                id1="Coherence",
                func_type=6,
                rsp_node=2,
                rsp_dir=3,
                ref_node=1,
                ref_dir=3,
                abscissa_spacing=1,
                abscissa_min=0.0,
                abscissa_inc=2.5,
                abscissa_spec_data_type=18,
                ordinate_spec_data_type=0,
                orddenom_spec_data_type=0,
                num_pts=len(freq),
                ord_data_type=2,
                x=freq,
                data=np.linspace(1.0, 0.5, len(freq)),
            ),
            mode="add",
        )
        print("wrote", path)

    uff, path = fresh("pyuff_frf_only.uff")
    uff.write_sets(frfs(0)[0], mode="add")
    print("wrote", path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--gmsh", action="store_true", help="also regenerate the gmsh files"
    )
    parser.add_argument(
        "--pyuff", action="store_true", help="also regenerate the pyuff files"
    )
    args = parser.parse_args()
    write_default()
    if args.gmsh:
        write_gmsh()
    if args.pyuff:
        write_pyuff()


if __name__ == "__main__":
    main()
