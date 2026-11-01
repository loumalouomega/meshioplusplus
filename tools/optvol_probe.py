#!/usr/bin/env python3
"""Measure ``optimize_volume`` on Netgen's tetrahedra (roadmap 3.2.1's probe).

Builds closed triangle skins with meshio++'s own primitives (an icosahedral
sphere and a capped cylinder), fills them with Netgen through its Python
package (``pip install netgen-mesher``; the ``netgen.meshing`` API, not the
nglib the ``tetrahedralize`` backend links) with and without Netgen's own
optimisation passes, and runs ``optimize_volume`` under both objectives. It
prints, per mesh, the least dihedral angle in and out and the tets below
15 degrees -- the table in ``doc/optimize_volume.md`` and
``doc/tetrahedralize.md``.

    python tools/optvol_probe.py [--maxh 0.2] [--json]

The probe holds when the ``worst`` column is never below the ``in`` column.
"""

from __future__ import annotations

import argparse
import json
import sys

import numpy as np

import meshioplusplus as mp


def skins():
    """The closed triangle skins: (name, mesh)."""
    sphere = mp.icosa_sphere(1.0, 3)
    solid = mp.cylinder(radius=1.0, height=2.0, nxy=12, nz=10, cell_type="wedge")
    cylinder = mp.convert_cells(mp.extract_skin(solid), mode="simplexify")
    return [("sphere", sphere), ("capped cylinder", cylinder)]


def netgen_volume(skin, maxh, optsteps):
    """Fill a closed triangle skin with Netgen's tetrahedra."""
    import netgen.meshing as nm

    mesh = nm.Mesh(dim=3)
    ids = [mesh.Add(nm.MeshPoint(nm.Pnt(*map(float, p)))) for p in skin.points]
    fd = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    for block in skin.cells:
        if block.type != "triangle":
            continue
        for tri in np.asarray(block.data):
            mesh.Add(nm.Element2D(fd, [ids[int(v)] for v in tri]))
    mesh.GenerateVolumeMesh(mp=nm.MeshingParameters(maxh=maxh, optsteps3d=optsteps))
    points = np.array([[p[0], p[1], p[2]] for p in mesh.Points()])
    tets = np.array([[v.nr - 1 for v in el.vertices] for el in mesh.Elements3D()])
    return mp.Mesh(points, [("tetra", tets)])


def dihedrals(mesh):
    q = mp.attach_quality(mesh).cell_data["quality:min_dihedral"]
    return np.concatenate([np.asarray(a, dtype=float).ravel() for a in q])


def measure(maxh=0.2):
    rows = []
    for name, skin in skins():
        for optsteps in (0, 3):
            vol = netgen_volume(skin, maxh, optsteps)
            row = {
                "skin": name,
                "netgen_passes": optsteps,
                "tets": len(vol.cells[0].data),
                "in": float(dihedrals(vol).min()),
                "slivers_in": int((dihedrals(vol) < 15.0).sum()),
            }
            for objective in ("energy", "worst"):
                out, rep = mp.optimize_volume(
                    vol, objective=objective, return_report=True
                )
                d = dihedrals(out)
                row[objective] = float(d.min())
                row["slivers_" + objective] = int((d < 15.0).sum())
                if objective == "worst":
                    row["moves"] = {
                        k: int(rep[k])
                        for k in (
                            "num_23_flips",
                            "num_32_flips",
                            "num_edge_removals",
                            "num_multi_face_removals",
                            "num_worst_smoothed",
                        )
                    }
            rows.append(row)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--maxh", type=float, default=0.2)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    try:
        import netgen.meshing  # noqa: F401
    except ImportError:
        print("optvol_probe: needs `pip install netgen-mesher`", file=sys.stderr)
        return 2
    rows = measure(args.maxh)
    if args.json:
        print(json.dumps(rows, indent=1))
    else:
        print(
            "| skin | Netgen passes | tets | least angle in | energy | worst "
            "| < 15 deg in / energy / worst |"
        )
        print("|---|---|---|---|---|---|---|")
        for r in rows:
            print(
                f"| {r['skin']} | {r['netgen_passes']} | {r['tets']} | {r['in']:.2f} "
                f"| {r['energy']:.2f} | {r['worst']:.2f} "
                f"| {r['slivers_in']} / {r['slivers_energy']} / {r['slivers_worst']} |"
            )
    return 0 if all(r["worst"] >= r["in"] for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
