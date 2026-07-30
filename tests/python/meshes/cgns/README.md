<!--pytest-codeblocks:skipfile-->
# CGNS reference fixtures

## `tri_quad_fields_cgnslib.cgns`

**Written by real cgnslib, not by meshio++.** This is the external-validation
layer `doc/formats/cgns.md` records as a follow-up of the v9.8.0 CGNS rewrite:
without a file some other implementation produced, CI can only ever prove that
meshio++ agrees with itself.

Provenance — reproduce with:

```bash
# cgnslib is on conda-forge (not in apt on a machine without sudo); this repo's
# other micromamba envs (xcpp, r-notebooks) live in the same gitignored root.
micromamba create -y -r .micromamba -n cgns -c conda-forge cgns   # 4.5.2
```

then, from the repo root with the package importable:

```python
import numpy as np, subprocess, tempfile
import meshioplusplus
from meshioplusplus import _core

pts = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0],
                [2, 1, 0], [1, 1, 0], [0, 1, 0]], float)
m = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 4], [0, 4, 5]])),
                              ("quad",     np.array([[1, 2, 3, 4]]))])
m.point_data["Temperature"] = np.arange(6, dtype=float) * 1.5
m.point_data["Velocity"] = np.arange(18, dtype=float).reshape(6, 3)
m.cell_data["Density"] = [np.array([1.0, 2.0]), np.array([3.0])]

ours = tempfile.mktemp(suffix=".cgns")
_core.cgns_write(ours, m, -1)
# cgnsconvert rewrites the file end to end through cgnslib's OWN API, so the
# committed bytes are cgnslib's, not ours -- that is the whole point.
subprocess.run([".micromamba/envs/cgns/bin/cgnsconvert", "-h", "-f", ours,
                "tests/python/meshes/cgns/tri_quad_fields_cgnslib.cgns"], check=True)
```

Produced with **cgnslib 4.5.2** (conda-forge `cgns-4.5.2-haeb83df_0`, linked
against HDF5 2.1.0). It exercises the whole v9.9.0 `FlowSolution_t` surface: two
cell blocks, a scalar and a 3-component vector `Vertex` field (so the
`Velocity_0/_1/_2` split-and-rejoin convention is covered), and a `CellCenter`
field spanning both blocks.

`cgnscheck -w3` reports **zero errors** on it, and on every mesh meshio++'s own
writer produces (see `test_cgnscheck_accepts_our_output`). It does report style
*warnings* — no `Family_t` on the zone, no `DataClass_t`, and "not a CGNS
data-name identifier" for any field whose name is not one of SIDS's standard
names (`Density`, `VelocityX`, ...). Those are recommendations, not
conformance failures, and meshio++ deliberately preserves the caller's own
field names rather than renaming them to fit SIDS's vocabulary.

Kept small (~18 KB) on purpose. Stored via Git LFS (`*.cgns` in
`.gitattributes`).
